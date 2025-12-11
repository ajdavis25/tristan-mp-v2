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

  ! phase 4A: uniform vertical seed field
  real    :: Omega           = 0.0
  real    :: q_shear         = 0.0
  real    :: backgr_T        = 0.0
  logical :: enable_rotation = .true.
  logical :: use_rotation    = .false.
  real    :: box_xcenter     = 0.0
  real    :: box_ycenter     = 0.0
  real    :: B0_ext          = 0.0

  private :: userSpatialDistribution
  private :: imprint_velocity_shear

contains

!=====================================================================
!  userReadInput
!=====================================================================
  subroutine userReadInput()
    implicit none

    call getInput('problem','Omega',           Omega,           0.0)
    call getInput('problem','q_shear',         q_shear,         0.0)
    call getInput('problem','backgr_T',        backgr_T,        0.0)
    call getInput('problem','enable_rotation', enable_rotation, .true.)
    call getInput('problem', 'B0_ext',         B0_ext,          0.0)

    ! phase 4A: enable full EM evolution now
    enable_fieldsolver    = .true.
    enable_currentdeposit = .true.

    box_xcenter = 0.5 * global_mesh%sx
#ifdef twoD
    box_ycenter = 0.5 * global_mesh%sy
#else
    box_ycenter = 0.0
#endif

    use_rotation = enable_rotation .and. (abs(Omega) > 0.0)

    print *, ">>> [3B] userReadInput:"
    print *, "    Omega          =", Omega
    print *, "    q_shear        =", q_shear
    print *, "    backgr_T       =", backgr_T
    print *, "    enable_rotation=", enable_rotation
    print *, "    box_xcenter    =", box_xcenter, " box_ycenter=", box_ycenter
    print *, "    use_rotation   =", use_rotation
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
! userInitParticles: thermal loading + optional shear imprint
!=====================================================================
  subroutine userInitParticles()
    implicit none
    integer :: s
    type(region) :: whole_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    spat_distr_ptr => userSpatialDistribution

    whole_region%x_min = 0.0
    whole_region%x_max = global_mesh%sx
#ifdef twoD
    whole_region%y_min = 0.0
    whole_region%y_max = global_mesh%sy
#endif

    print *, ">>> [3B] userInitParticles: thermal loading (zero current)"

    do s = 1, size(species)
      call fillRegionWithThermalPlasma(whole_region, (/s/), 1, ppc0, backgr_T, &
                                       zero_current=.true., spat_distr_ptr=spat_distr_ptr)
    end do

    if (use_rotation .and. (q_shear /= 0.0)) then
      print *, ">>> [3B] userInitParticles: imprinting shearing flow:"
      print *, "    v_y(x) = -q_shear * Omega * (x - x_center)"
      call imprint_velocity_shear()
    else
      print *, ">>> [3B] userInitParticles: rotation/shear disabled (no v_y profile)"
    end if
  end subroutine userInitParticles

!=====================================================================
! imprint_velocity_shear: v_y(x) = -q_shear * Omega * (x - box_xcenter)
!=====================================================================
  subroutine imprint_velocity_shear()
    implicit none
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i_idx
    real :: x_glob
    real :: u, v, w
    real :: vx, vy_new, vz
    real :: gamma_old, gamma_new, beta2, scale
    real :: u_new, v_new, w_new
    real :: vy_shear

    print *, ">>> [3B] imprint_velocity_shear: starting shear loop..."

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz

        do p = 1, species(s)%prtl_tile(ti,tj,tk)%npart_sp

          ! cell index -> global X
          i_idx  = species(s)%prtl_tile(ti,tj,tk)%xi(p)
          x_glob = real(i_idx + this_meshblock%ptr%x0) + &
                   species(s)%prtl_tile(ti,tj,tk)%dx(p)

          ! target shearing-sheet drift:
          !   v_y,sh(x) = - q_shear * Omega * (x - box_xcenter)
          vy_shear = -q_shear * Omega * (x_glob - box_xcenter)

          ! old momenta
          u = species(s)%prtl_tile(ti,tj,tk)%u(p)
          v = species(s)%prtl_tile(ti,tj,tk)%v(p)
          w = species(s)%prtl_tile(ti,tj,tk)%w(p)

          ! convert to velocities
          gamma_old = sqrt(1.0 + u*u + v*v + w*w)
          ! freeze x-motion to avoid mixing away the profile
          vx = 0.0
          vz = CC * w / gamma_old

          ! overwrite v_y with shearing-sheet drift
          vy_new = vy_shear

          ! beta^2 safety check
          beta2 = (vx*vx + vy_new*vy_new + vz*vz) / (CC*CC)
          if (beta2 >= 1.0) then
            scale  = 0.9 / sqrt(beta2)
            vx     = vx     * scale
            vy_new = vy_new * scale
            vz     = vz     * scale
            beta2  = (vx*vx + vy_new*vy_new + vz*vz) / (CC*CC)
          end if

          ! back to momenta
          gamma_new = 1.0 / sqrt(1.0 - beta2)
          u_new = gamma_new * vx     / CC
          v_new = gamma_new * vy_new / CC
          w_new = gamma_new * vz     / CC

          species(s)%prtl_tile(ti,tj,tk)%u(p) = u_new
          species(s)%prtl_tile(ti,tj,tk)%v(p) = v_new
          species(s)%prtl_tile(ti,tj,tk)%w(p) = w_new

        end do

      end do
      end do
      end do
    end do

    print *, ">>> [3B] imprint_velocity_shear: DONE."
  end subroutine imprint_velocity_shear

!=====================================================================
! userInitFields: zero fields
!=====================================================================
  subroutine userInitFields()
    implicit none
    ex(:,:,:) = 0.0; ey(:,:,:) = 0.0; ez(:,:,:) = 0.0

    bx(:,:,:) = 0.0
    by(:,:,:) = 0.0
    bz(:,:,:) = B0_ext  ! uniform vertical field

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
! userDriveParticles: no extra body forces (just lorentz from Bz)
!--------------------------------------------------------------------
  subroutine userDriveParticles(step)
    implicit none
    integer, intent(in), optional :: step
    ! no external forces – lorentz force comes from main pusher
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
