module m_userfile
#if defined(MPI08) || defined(MPI)
  use mpi
#endif
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
#ifdef USROUTPUT
  use m_writeusroutput
#endif
  implicit none

  ! ---------- define PI constant ----------
  real, parameter :: PI_F = acos(-1.0)
  ! ----------------------------------------

  ! phase 4A/4B: shear + uniform vertical seed field + optional perturbation
  real    :: Omega             = 0.0
  real    :: q_shear           = 0.0
  real    :: backgr_T          = 0.0
  logical :: enable_rotation   = .true.
  logical :: use_rotation      = .false.
  real    :: box_xcenter       = 0.0
  real    :: box_ycenter       = 0.0
  real    :: B0_ext            = 0.0

  ! phase 4B: small velocity perturbation to seed MRI-like modes
  logical :: enable_perturbation = .false.
  real    :: eps_pert            = 0.0   ! dimensionless amplitude (fraction of c)
  integer :: ntimes_user         = 1     ! particle sub-cycles per global timestep
  integer, save :: drive_last_step = -HUGE(1)
  logical, save :: drive_dt_mode_set = .false.
  logical, save :: drive_single_call = .false.

  private :: userSpatialDistribution
  private :: imprint_velocity_shear
  private :: imprint_perturbation

contains

!=====================================================================
!  userReadInput
!=====================================================================
  subroutine userReadInput()
    implicit none

    call getInput('problem','Omega',              Omega,              0.0)
    call getInput('problem','q_shear',            q_shear,            0.0)
    call getInput('problem','backgr_T',           backgr_T,           0.0)
    call getInput('problem','enable_rotation',    enable_rotation,    .true.)
    call getInput('problem','B0_ext',             B0_ext,             0.0)
    call getInput('problem','enable_perturbation',enable_perturbation,.false.)
    call getInput('problem','eps_pert',           eps_pert,           0.0)
    call getInput('algorithm','ntimes',           ntimes_user,         1)

    box_xcenter = 0.5 * global_mesh%sx
#ifdef twoD
    box_ycenter = 0.5 * global_mesh%sy
#else
    box_ycenter = 0.0
#endif

    use_rotation = enable_rotation .and. (abs(Omega) > 0.0)

    print *, ">>> [4] userReadInput:"
    print *, "    Omega              =", Omega
    print *, "    q_shear            =", q_shear
    print *, "    backgr_T           =", backgr_T
    print *, "    enable_rotation    =", enable_rotation
    print *, "    B0_ext             =", B0_ext
    print *, "    enable_perturbation=", enable_perturbation
    print *, "    eps_pert           =", eps_pert
    print *, "    ntimes_user        =", ntimes_user
    print *, "    box_xcenter        =", box_xcenter, " box_ycenter=", box_ycenter
    print *, "    use_rotation       =", use_rotation
  end subroutine userReadInput

!=====================================================================
! spatial PDF
!=====================================================================
  function userSpatialDistribution(x_glob,y_glob,z_glob,d1,d2,d3) result(weight)
    real :: weight
    real, intent(in), optional :: x_glob,y_glob,z_glob,d1,d2,d3
    weight = 1.0
  end function userSpatialDistribution


!=====================================================================
! userInitParticles
!=====================================================================
subroutine userInitParticles()
  implicit none
  integer :: s
  type(region) :: whole_region
  procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

  spat_distr_ptr => userSpatialDistribution

  whole_region%x_min = 0.0
  whole_region%x_max = global_mesh%sx
#if defined(twoD) || defined(threeD)
  whole_region%y_min = 0.0
  whole_region%y_max = global_mesh%sy
#endif
#if defined(threeD)
  whole_region%z_min = 0.0
  whole_region%z_max = global_mesh%sz
#endif

  print *, ">>> [4] userInitParticles: thermal loading (zero current)"

  do s = 1, size(species)
    call fillRegionWithThermalPlasma(whole_region, (/s/), 1, ppc0, backgr_T, &
                                     zero_current=.true., spat_distr_ptr=spat_distr_ptr)
  end do

  ! For MRI runs, the background Keplerian shear is the equilibrium that
  ! cancels the tidal term in Hill's equations. Without this, the tidal
  ! forcing accelerates the whole box and can drive large, unphysical
  ! transients before the instability is even seeded.
  if (use_rotation .and. (q_shear /= 0.0)) then
    print *, ">>> [4] imprinting shearing flow"
    call imprint_velocity_shear()
  else
    print *, ">>> [4] rotation/shear disabled"
  end if

  if (enable_perturbation .and. (eps_pert /= 0.0)) then
    print *, ">>> [4B] imprinting dvx perturbation"
    call imprint_perturbation()
  else
    print *, ">>> [4B] no perturbation applied"
  end if
end subroutine userInitParticles


!=====================================================================
! imprint_velocity_shear
!=====================================================================
  subroutine imprint_velocity_shear()
    implicit none
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i_idx
    real :: x_glob, u, v, w
    real :: vx, vy_new, vz
    real :: gamma_old, gamma_new, beta2, scale
    real :: u_new, v_new, w_new
    real :: vy_shear

    print *, ">>> [3B] imprint_velocity_shear..."

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti,tj,tk)%npart_sp

          i_idx  = species(s)%prtl_tile(ti,tj,tk)%xi(p)
          x_glob = real(i_idx + this_meshblock%ptr%x0) + species(s)%prtl_tile(ti,tj,tk)%dx(p)

          vy_shear = -q_shear * Omega * (x_glob - box_xcenter)

          u = species(s)%prtl_tile(ti,tj,tk)%u(p)
          v = species(s)%prtl_tile(ti,tj,tk)%v(p)
          w = species(s)%prtl_tile(ti,tj,tk)%w(p)

          gamma_old = sqrt(1.0 + u*u + v*v + w*w)
          vx = 0.0
          vz = CC * w / gamma_old
          vy_new = vy_shear

          beta2 = (vx*vx + vy_new*vy_new + vz*vz)/(CC*CC)
          if (beta2 >= 1.0) then
            scale = 0.9/sqrt(beta2)
            vx = vx*scale; vy_new = vy_new*scale; vz = vz*scale
            beta2 = (vx*vx + vy_new*vy_new + vz*vz)/(CC*CC)
          end if

          gamma_new = 1.0/sqrt(1.0 - beta2)
          u_new = gamma_new * vx / CC
          v_new = gamma_new * vy_new / CC
          w_new = gamma_new * vz / CC

          species(s)%prtl_tile(ti,tj,tk)%u(p) = u_new
          species(s)%prtl_tile(ti,tj,tk)%v(p) = v_new
          species(s)%prtl_tile(ti,tj,tk)%w(p) = w_new

        end do
      end do
      end do
      end do
    end do

    print *, ">>> [3B] imprint_velocity_shear DONE"
  end subroutine imprint_velocity_shear

!=====================================================================
! imprint_perturbation  (MRI seed mode)
!=====================================================================
  subroutine imprint_perturbation()
    implicit none
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: k_idx
    real :: z_glob, u, v, w
    real :: vx, vy, vz, dvx, dvy
    real :: gamma_old, gamma_new, beta2, scale
    real :: u_new, v_new, w_new
    real :: kz
    logical :: perturb_com

    if (.not. enable_perturbation) return

    print *, ">>> [4B] imprint_perturbation..."

    ! Epicyclic validation: excite a pure center-of-mass kick (no k-modes).
    ! Keep the original k=1 sinusoidal seed for MRI-like runs with B0_ext /= 0.
    if (abs(B0_ext) > 0.0) then
      kz = 2.0 * PI_F / global_mesh%sz
      perturb_com = .false.
      print *, ">>> [4B] perturbation: k=1 sin(kz*z), cos(kz*z) seed"
    else
      kz = 0.0
      perturb_com = .true.
      print *, ">>> [4B] perturbation: COM dvx kick"
    end if

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti,tj,tk)%npart_sp

          k_idx  = species(s)%prtl_tile(ti,tj,tk)%zi(p)
          z_glob = real(k_idx + this_meshblock%ptr%z0) + species(s)%prtl_tile(ti,tj,tk)%dz(p)

          u = species(s)%prtl_tile(ti,tj,tk)%u(p)
          v = species(s)%prtl_tile(ti,tj,tk)%v(p)
          w = species(s)%prtl_tile(ti,tj,tk)%w(p)

          gamma_old = sqrt(1.0 + u*u + v*v + w*w)
          vx = CC * u / gamma_old
          vy = CC * v / gamma_old
          vz = CC * w / gamma_old

          if (perturb_com) then
            dvx = eps_pert * CC
            dvy = 0.0
          else
            dvx = eps_pert * CC * sin(kz * z_glob)
            dvy = eps_pert * CC * cos(kz * z_glob)
          end if
          vx  = vx + dvx
          vy  = vy + dvy

          beta2 = (vx*vx + vy*vy + vz*vz)/(CC*CC)
          if (beta2 >= 1.0) then
            scale = 0.9/sqrt(beta2)
            vx = vx*scale; vy = vy*scale; vz = vz*scale
            beta2 = (vx*vx + vy*vy + vz*vz)/(CC*CC)
          end if

          gamma_new = 1.0/sqrt(1.0 - beta2)
          u_new = gamma_new * vx / CC
          v_new = gamma_new * vy / CC
          w_new = gamma_new * vz / CC

          species(s)%prtl_tile(ti,tj,tk)%u(p) = u_new
          species(s)%prtl_tile(ti,tj,tk)%v(p) = v_new
          species(s)%prtl_tile(ti,tj,tk)%w(p) = w_new

        end do
      end do
      end do
      end do
    end do

    print *, ">>> [4B] imprint_perturbation DONE"
  end subroutine imprint_perturbation

!=====================================================================
! userInitFields
!=====================================================================
  subroutine userInitFields()
    implicit none
    ex(:,:,:) = 0.0; ey(:,:,:) = 0.0; ez(:,:,:) = 0.0
    bx(:,:,:) = 0.0; by(:,:,:) = 0.0
    bz(:,:,:) = B0_ext
    jx(:,:,:) = 0.0; jy(:,:,:) = 0.0; jz(:,:,:) = 0.0
    print *, ">>> [4A] userInitFields initialized Bz =", B0_ext
  end subroutine userInitFields

!--------------------------------------------------------------------
! userCurrentDeposit: required by tristanmainloop.F90
!--------------------------------------------------------------------
  subroutine userCurrentDeposit(step)
    implicit none
    integer, intent(in), optional :: step
    ! no custom current sources
  end subroutine userCurrentDeposit

!--------------------------------------------------------------------
! userDriveParticles: Coriolis + tidal forces (Hill/shearing-box)
!
! notes:
!  - `userDriveParticles()` is called exactly once per global timestep from
!    `src/main/tristanmainloop.F90` (not inside the mover substeps).
!  - in TRISTAN code units, the timestep is dt=1; CC sets c*dt/dx.
!  - if this routine modifies momenta after `moveParticles()`, it must also
!    adjust particle positions so charge-conserving current deposition remains
!    consistent (see below).
!--------------------------------------------------------------------
subroutine userDriveParticles(step)
  implicit none
  integer, intent(in), optional :: step

  integer :: s, ti, tj, tk, p
  integer(kind=2) :: i_idx
  real :: x_glob, x_rel
  real :: u, v, w
  real :: vx, vy, vz
  real :: vx_old, vy_old, vz_old
  real :: gamma_old, gamma_new
  real :: ax, ay
  real :: v2, scale
  real :: dt_loc
  integer(kind=2) :: temp_i
  real :: temp_r
  real, parameter :: beta2_max = 1.0 - 1.0e-6

  if (.not. use_rotation) return

  ! in TRISTAN code units, `step` advances by 1 each global timestep.
  dt_loc = 1.0

  do s = 1, nspec
    do ti = 1, species(s)%tile_nx
    do tj = 1, species(s)%tile_ny
    do tk = 1, species(s)%tile_nz
      do p = 1, species(s)%prtl_tile(ti,tj,tk)%npart_sp

        ! --- particle position (cell units) ---
        i_idx  = species(s)%prtl_tile(ti,tj,tk)%xi(p)
        x_glob = real(i_idx + this_meshblock%ptr%x0) + &
                 species(s)%prtl_tile(ti,tj,tk)%dx(p)
        x_rel = x_glob - box_xcenter

        ! --- momenta ---
        u = species(s)%prtl_tile(ti,tj,tk)%u(p)
        v = species(s)%prtl_tile(ti,tj,tk)%v(p)
        w = species(s)%prtl_tile(ti,tj,tk)%w(p)

        gamma_old = sqrt(1.0 + u*u + v*v + w*w)

        ! --- velocities ---
        vx_old = CC * u / gamma_old
        vy_old = CC * v / gamma_old
        vz_old = CC * w / gamma_old
        vx = vx_old
        vy = vy_old
        vz = vz_old

        !------------------------------------------------------------
        ! hill's equations (shearing box, rotating frame):
        !   dvx/dt =  2Ω vy + 2 q Ω^2 (x - x_center)
        !   dvy/dt = -2Ω vx
        !   dvz/dt =  0
        !------------------------------------------------------------
        ax =  2.0 * Omega * vy + 2.0 * q_shear * Omega * Omega * x_rel
        ay = -2.0 * Omega * vx

        ! --- update velocities (explicit) ---
        vx = vx + ax * dt_loc
        vy = vy + ay * dt_loc
        ! vz unchanged

        ! --- enforce |v| < c ---
        v2 = vx*vx + vy*vy + vz*vz
        if (v2 >= beta2_max * CC * CC) then
          scale = sqrt((beta2_max * CC * CC) / v2)
          vx = vx * scale
          vy = vy * scale
          vz = vz * scale
          v2 = vx*vx + vy*vy + vz*vz
        end if

        ! --- back to momentum ---
        gamma_new = 1.0 / sqrt(1.0 - v2/(CC*CC))
        species(s)%prtl_tile(ti,tj,tk)%u(p) = gamma_new * vx / CC
        species(s)%prtl_tile(ti,tj,tk)%v(p) = gamma_new * vy / CC
        species(s)%prtl_tile(ti,tj,tk)%w(p) = gamma_new * vz / CC

        !------------------------------------------------------------
        ! IMPORTANT for charge-conserving current deposition:
        ! `depositCurrents()` reconstructs the start-of-step position by
        !   x1 = x2 - v_end * dt
        ! so if we modify momenta after `moveParticles()`, we must also
        ! adjust the particle position so that the implied displacement
        ! uses the same v_end. Otherwise we inject nonphysical currents.
        !
        ! here dt=1 (code units), so shift the *final* position by
        !   Δx = (v_new - v_old) * dt
        ! using the same carry logic as `position_update.F08`.
        !------------------------------------------------------------
#if defined(oneD) || defined(twoD) || defined(threeD)
        species(s)%prtl_tile(ti,tj,tk)%dx(p) = species(s)%prtl_tile(ti,tj,tk)%dx(p) + (vx - vx_old) * dt_loc
        temp_i = INT(species(s)%prtl_tile(ti,tj,tk)%dx(p), 2)
        temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti,tj,tk)%dx(p)) + temp_i, REAL(temp_i)) - 1
        temp_i = INT(temp_r, 2)
        species(s)%prtl_tile(ti,tj,tk)%xi(p) = species(s)%prtl_tile(ti,tj,tk)%xi(p) + temp_i
        species(s)%prtl_tile(ti,tj,tk)%dx(p) = species(s)%prtl_tile(ti,tj,tk)%dx(p) - temp_r
#endif

#if defined(twoD) || defined(threeD)
        species(s)%prtl_tile(ti,tj,tk)%dy(p) = species(s)%prtl_tile(ti,tj,tk)%dy(p) + (vy - vy_old) * dt_loc
        temp_i = INT(species(s)%prtl_tile(ti,tj,tk)%dy(p), 2)
        temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti,tj,tk)%dy(p)) + temp_i, REAL(temp_i)) - 1
        temp_i = INT(temp_r, 2)
        species(s)%prtl_tile(ti,tj,tk)%yi(p) = species(s)%prtl_tile(ti,tj,tk)%yi(p) + temp_i
        species(s)%prtl_tile(ti,tj,tk)%dy(p) = species(s)%prtl_tile(ti,tj,tk)%dy(p) - temp_r
#endif

#if defined(threeD)
        species(s)%prtl_tile(ti,tj,tk)%dz(p) = species(s)%prtl_tile(ti,tj,tk)%dz(p) + (vz - vz_old) * dt_loc
        temp_i = INT(species(s)%prtl_tile(ti,tj,tk)%dz(p), 2)
        temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti,tj,tk)%dz(p)) + temp_i, REAL(temp_i)) - 1
        temp_i = INT(temp_r, 2)
        species(s)%prtl_tile(ti,tj,tk)%zi(p) = species(s)%prtl_tile(ti,tj,tk)%zi(p) + temp_i
        species(s)%prtl_tile(ti,tj,tk)%dz(p) = species(s)%prtl_tile(ti,tj,tk)%dz(p) - temp_r
#endif

      end do
    end do
    end do
    end do
  end do

end subroutine userDriveParticles


!--------------------------------------------------------------------
! userFieldBoundaryConditions: required by tristanmainloop.F90
!--------------------------------------------------------------------
  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, intent(in), optional :: step
    logical, intent(in), optional :: updateE, updateB
    ! no field BC modifications
  end subroutine userFieldBoundaryConditions

!--------------------------------------------------------------------
! userParticleBoundaryConditions: required by tristanmainloop.F90
!--------------------------------------------------------------------
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, intent(in), optional :: step
    ! no particle BC modifications
  end subroutine userParticleBoundaryConditions

!--------------------------------------------------------------------
! userDeallocate: required by finalize.F90
!--------------------------------------------------------------------
  subroutine userDeallocate()
    implicit none
    ! nothing to free
  end subroutine userDeallocate

!--------------------------------------------------------------------
! writeUsrRestart / readUsrRestart: required by restart.F90
!--------------------------------------------------------------------
  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
    ! nothing to write
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
    ! nothing to read
  end subroutine readUsrRestart

end module m_userfile
