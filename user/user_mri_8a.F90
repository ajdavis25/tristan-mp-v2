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

  real, parameter :: PI_F = acos(-1.0)
  integer, parameter :: N_TRANSITION_SLABS = 4

  real    :: Omega                 = 0.0
  real    :: q_shear               = 0.0
  real    :: B0_ext                = 0.0
  real    :: backgr_T_inner        = 1.0
  real    :: backgr_T_outer        = 1.0
  real    :: inner_xmin_frac       = 0.40
  real    :: inner_xmax_frac       = 0.60
  real    :: transition_width_frac = 0.05
  real    :: buffer_width_frac     = 0.125
  real    :: field_damp_kappa      = 10.0
  real    :: particle_damp_kappa   = 8.0
  real    :: particle_kill_width_frac = 0.03125
  logical :: enable_rotation       = .true.
  logical :: use_rotation          = .false.
  integer :: drive_inner_only      = 1
  logical :: enable_perturbation   = .false.
  integer :: perturb_inner_only    = 1
  real    :: eps_pert              = 0.0
  integer :: ntimes_user           = 1

  real    :: box_xcenter           = 0.0
  real    :: box_ycenter           = 0.0
  real    :: sx_glob               = 0.0
  real    :: zone_x0               = 0.0
  real    :: zone_x1               = 0.0
  real    :: zone_x2               = 0.0
  real    :: zone_x3               = 0.0
  real    :: buffer_edge           = 0.0
  real    :: kill_edge             = 0.0
  logical :: boundary_reported     = .false.

  private :: userSpatialDistribution
  private :: imprint_velocity_shear
  private :: imprint_perturbation
  private :: validate_zone_geometry
  private :: build_zone_edges
  private :: load_temperature_slab
  private :: load_transition_slabs
  private :: x_drive_window
  private :: x_seed_window
  private :: x_buffer_depth

contains

  subroutine userReadInput()
    implicit none

    call getInput('problem', 'Omega',                    Omega,                    0.0)
    call getInput('problem', 'q_shear',                  q_shear,                  0.0)
    call getInput('problem', 'B0_ext',                   B0_ext,                   0.0)
    call getInput('problem', 'backgr_T_inner',           backgr_T_inner,           1.0)
    call getInput('problem', 'backgr_T_outer',           backgr_T_outer,           1.0)
    call getInput('problem', 'inner_xmin_frac',          inner_xmin_frac,          0.40)
    call getInput('problem', 'inner_xmax_frac',          inner_xmax_frac,          0.60)
    call getInput('problem', 'transition_width_frac',    transition_width_frac,    0.05)
    call getInput('problem', 'buffer_width_frac',        buffer_width_frac,        0.125)
    call getInput('problem', 'field_damp_kappa',         field_damp_kappa,         10.0)
    call getInput('problem', 'particle_damp_kappa',      particle_damp_kappa,      8.0)
    call getInput('problem', 'particle_kill_width_frac', particle_kill_width_frac, 0.03125)
    call getInput('problem', 'enable_rotation',          enable_rotation,          .true.)
    call getInput('problem', 'drive_inner_only',         drive_inner_only,         1)
    call getInput('problem', 'enable_perturbation',      enable_perturbation,      .false.)
    call getInput('problem', 'perturb_inner_only',       perturb_inner_only,       1)
    call getInput('problem', 'eps_pert',                 eps_pert,                 0.0)
    call getInput('algorithm', 'ntimes',                 ntimes_user,              1)

    sx_glob = real(global_mesh%sx)
    box_xcenter = 0.5 * sx_glob
#ifdef twoD
    box_ycenter = 0.5 * real(global_mesh%sy)
#else
    box_ycenter = 0.0
#endif

    use_rotation = enable_rotation .and. (abs(Omega) > 0.0)

    call validate_zone_geometry()
    call build_zone_edges()

    if (mpi_rank .eq. 0) then
      print *, ">>> [8A] userReadInput:"
      print *, "    Omega                    =", Omega
      print *, "    q_shear                  =", q_shear
      print *, "    B0_ext                   =", B0_ext
      print *, "    backgr_T_inner           =", backgr_T_inner
      print *, "    backgr_T_outer           =", backgr_T_outer
      print *, "    inner_xmin_frac          =", inner_xmin_frac
      print *, "    inner_xmax_frac          =", inner_xmax_frac
      print *, "    transition_width_frac    =", transition_width_frac
      print *, "    buffer_width_frac        =", buffer_width_frac
      print *, "    field_damp_kappa         =", field_damp_kappa
      print *, "    particle_damp_kappa      =", particle_damp_kappa
      print *, "    particle_kill_width_frac =", particle_kill_width_frac
      print *, "    enable_rotation          =", enable_rotation
      print *, "    drive_inner_only         =", drive_inner_only
      print *, "    enable_perturbation      =", enable_perturbation
      print *, "    perturb_inner_only       =", perturb_inner_only
      print *, "    eps_pert                 =", eps_pert
      print *, "    ntimes_user              =", ntimes_user
      print *, "    box_xcenter              =", box_xcenter, " box_ycenter=", box_ycenter
      print *, "    use_rotation             =", use_rotation
      print *, "    zone_x0                  =", zone_x0
      print *, "    zone_x1                  =", zone_x1
      print *, "    zone_x2                  =", zone_x2
      print *, "    zone_x3                  =", zone_x3
      print *, "    buffer_edge              =", buffer_edge
      print *, "    kill_edge                =", kill_edge
    end if
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, d1, d2, d3) result(weight)
    real :: weight
    real, intent(in), optional :: x_glob, y_glob, z_glob, d1, d2, d3
    weight = 1.0
  end function userSpatialDistribution

  subroutine validate_zone_geometry()
    implicit none
    real :: zone_x0_frac, zone_x3_frac

    if (inner_xmin_frac <= 0.0) then
      print *, "ERROR [8A]: inner_xmin_frac must be > 0"
      stop
    end if
    if (inner_xmax_frac >= 1.0) then
      print *, "ERROR [8A]: inner_xmax_frac must be < 1"
      stop
    end if
    if (inner_xmin_frac >= inner_xmax_frac) then
      print *, "ERROR [8A]: inner_xmin_frac must be < inner_xmax_frac"
      stop
    end if
    if (transition_width_frac <= 0.0) then
      print *, "ERROR [8A]: transition_width_frac must be > 0"
      stop
    end if
    if (buffer_width_frac <= 0.0) then
      print *, "ERROR [8A]: buffer_width_frac must be > 0"
      stop
    end if
    if (particle_kill_width_frac <= 0.0) then
      print *, "ERROR [8A]: particle_kill_width_frac must be > 0"
      stop
    end if
    if (particle_kill_width_frac >= buffer_width_frac) then
      print *, "ERROR [8A]: particle_kill_width_frac must be < buffer_width_frac"
      stop
    end if

    zone_x0_frac = inner_xmin_frac - transition_width_frac
    zone_x3_frac = inner_xmax_frac + transition_width_frac
    if (zone_x0_frac < 0.0) then
      print *, "ERROR [8A]: left transition extends outside the box"
      stop
    end if
    if (zone_x3_frac > 1.0) then
      print *, "ERROR [8A]: right transition extends outside the box"
      stop
    end if
    if (buffer_width_frac >= zone_x0_frac) then
      print *, "ERROR [8A]: left buffer overlaps the source transition region"
      stop
    end if
    if ((1.0 - buffer_width_frac) <= zone_x3_frac) then
      print *, "ERROR [8A]: right buffer overlaps the source transition region"
      stop
    end if
  end subroutine validate_zone_geometry

  subroutine build_zone_edges()
    implicit none

    zone_x0 = sx_glob * (inner_xmin_frac - transition_width_frac)
    zone_x1 = sx_glob * inner_xmin_frac
    zone_x2 = sx_glob * inner_xmax_frac
    zone_x3 = sx_glob * (inner_xmax_frac + transition_width_frac)
    buffer_edge = sx_glob * buffer_width_frac
    kill_edge = sx_glob * particle_kill_width_frac
  end subroutine build_zone_edges

  subroutine load_temperature_slab(x_min, x_max, slab_temp)
    implicit none
    real, intent(in) :: x_min, x_max, slab_temp
    integer :: s
    type(region) :: fill_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    if (x_max <= x_min) return

    spat_distr_ptr => userSpatialDistribution

    fill_region%x_min = x_min
    fill_region%x_max = x_max
#if defined(twoD) || defined(threeD)
    fill_region%y_min = 0.0
    fill_region%y_max = real(global_mesh%sy)
#endif
#if defined(threeD)
    fill_region%z_min = 0.0
    fill_region%z_max = real(global_mesh%sz)
#endif

    do s = 1, nspec
      call fillRegionWithThermalPlasma(fill_region, (/s/), 1, ppc0, slab_temp, &
                                       zero_current=.true., spat_distr_ptr=spat_distr_ptr)
    end do
  end subroutine load_temperature_slab

  subroutine load_transition_slabs(x_left, x_right, temp_left, temp_right, label)
    implicit none
    real, intent(in) :: x_left, x_right, temp_left, temp_right
    character(len=*), intent(in) :: label
    integer :: i
    real :: slab_width, x_lo, x_hi, x_mid, alpha, slab_temp

    if (x_right <= x_left) return

    slab_width = (x_right - x_left) / real(N_TRANSITION_SLABS)
    do i = 1, N_TRANSITION_SLABS
      x_lo = x_left + real(i - 1) * slab_width
      x_hi = x_left + real(i) * slab_width
      x_mid = 0.5 * (x_lo + x_hi)
      alpha = (x_mid - x_left) / (x_right - x_left)
      slab_temp = temp_left + alpha * (temp_right - temp_left)
      if (mpi_rank .eq. 0) then
        print *, ">>> [8A] ", trim(label), " slab", i, " x=[", x_lo, ",", x_hi, "] T=", slab_temp
      end if
      call load_temperature_slab(x_lo, x_hi, slab_temp)
    end do
  end subroutine load_transition_slabs

  function x_drive_window(x_glob) result(window)
    implicit none
    real, intent(in) :: x_glob
    real :: window
    real :: arg

    if (drive_inner_only == 0) then
      window = 1.0
      return
    end if

    if ((x_glob < zone_x0) .or. (x_glob > zone_x3)) then
      window = 0.0
    else if (x_glob < zone_x1) then
      arg = PI_F * (x_glob - zone_x0) / (zone_x1 - zone_x0)
      window = 0.5 * (1.0 - cos(arg))
    else if (x_glob <= zone_x2) then
      window = 1.0
    else if (x_glob <= zone_x3) then
      arg = PI_F * (x_glob - zone_x2) / (zone_x3 - zone_x2)
      window = 0.5 * (1.0 + cos(arg))
    else
      window = 0.0
    end if
  end function x_drive_window

  function x_seed_window(x_glob) result(window)
    implicit none
    real, intent(in) :: x_glob
    real :: window
    real :: arg

    if (perturb_inner_only == 0) then
      window = 1.0
      return
    end if

    if ((x_glob < zone_x0) .or. (x_glob > zone_x3)) then
      window = 0.0
    else if (x_glob < zone_x1) then
      arg = PI_F * (x_glob - zone_x0) / (zone_x1 - zone_x0)
      window = 0.5 * (1.0 - cos(arg))
    else if (x_glob <= zone_x2) then
      window = 1.0
    else if (x_glob <= zone_x3) then
      arg = PI_F * (x_glob - zone_x2) / (zone_x3 - zone_x2)
      window = 0.5 * (1.0 + cos(arg))
    else
      window = 0.0
    end if
  end function x_seed_window

  function x_buffer_depth(x_glob) result(eta)
    implicit none
    real, intent(in) :: x_glob
    real :: eta
    real :: right_start

    right_start = sx_glob - buffer_edge
    eta = 0.0

    if (buffer_edge <= 0.0) return

    if (x_glob <= buffer_edge) then
      eta = (buffer_edge - x_glob) / buffer_edge
    else if (x_glob >= right_start) then
      eta = (x_glob - right_start) / buffer_edge
    end if

    eta = min(1.0, max(0.0, eta))
  end function x_buffer_depth

  subroutine userInitParticles()
    implicit none

    if (mpi_rank .eq. 0) then
      print *, ">>> [8A] userInitParticles: loading source-plus-buffer box"
      print *, ">>> [8A] outer-left core x=[0,", zone_x0, "] T=", backgr_T_outer
    end if
    call load_temperature_slab(0.0, zone_x0, backgr_T_outer)

    call load_transition_slabs(zone_x0, zone_x1, backgr_T_outer, backgr_T_inner, "left transition")

    if (mpi_rank .eq. 0) then
      print *, ">>> [8A] inner core x=[", zone_x1, ",", zone_x2, "] T=", backgr_T_inner
    end if
    call load_temperature_slab(zone_x1, zone_x2, backgr_T_inner)

    call load_transition_slabs(zone_x2, zone_x3, backgr_T_inner, backgr_T_outer, "right transition")

    if (mpi_rank .eq. 0) then
      print *, ">>> [8A] outer-right core x=[", zone_x3, ",", sx_glob, "] T=", backgr_T_outer
    end if
    call load_temperature_slab(zone_x3, sx_glob, backgr_T_outer)

    if (use_rotation .and. (q_shear /= 0.0)) then
      if (mpi_rank .eq. 0) print *, ">>> [8A] imprinting localized shearing flow"
      call imprint_velocity_shear()
    else
      if (mpi_rank .eq. 0) print *, ">>> [8A] rotation/shear disabled"
    end if

    if (enable_perturbation .and. (eps_pert /= 0.0)) then
      if (mpi_rank .eq. 0) print *, ">>> [8A] imprinting seeded perturbation"
      call imprint_perturbation()
    else
      if (mpi_rank .eq. 0) print *, ">>> [8A] no perturbation applied"
    end if
  end subroutine userInitParticles

  subroutine imprint_velocity_shear()
    implicit none
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i_idx
    real :: x_glob, u, v, w
    real :: vx, vy_new, vz
    real :: gamma_old, gamma_new, beta2, scale
    real :: u_new, v_new, w_new
    real :: vy_shear, drive_window

    if (mpi_rank .eq. 0) print *, ">>> [8A] imprint_velocity_shear..."

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp

          i_idx = species(s)%prtl_tile(ti, tj, tk)%xi(p)
          x_glob = real(i_idx + this_meshblock%ptr%x0) + species(s)%prtl_tile(ti, tj, tk)%dx(p)

          drive_window = x_drive_window(x_glob)
          vy_shear = -q_shear * Omega * (x_glob - box_xcenter) * drive_window

          u = species(s)%prtl_tile(ti, tj, tk)%u(p)
          v = species(s)%prtl_tile(ti, tj, tk)%v(p)
          w = species(s)%prtl_tile(ti, tj, tk)%w(p)

          gamma_old = sqrt(1.0 + u * u + v * v + w * w)
          vx = 0.0
          vz = CC * w / gamma_old
          vy_new = vy_shear

          beta2 = (vx * vx + vy_new * vy_new + vz * vz) / (CC * CC)
          if (beta2 >= 1.0) then
            scale = 0.9 / sqrt(beta2)
            vx = vx * scale
            vy_new = vy_new * scale
            vz = vz * scale
            beta2 = (vx * vx + vy_new * vy_new + vz * vz) / (CC * CC)
          end if

          gamma_new = 1.0 / sqrt(1.0 - beta2)
          u_new = gamma_new * vx / CC
          v_new = gamma_new * vy_new / CC
          w_new = gamma_new * vz / CC

          species(s)%prtl_tile(ti, tj, tk)%u(p) = u_new
          species(s)%prtl_tile(ti, tj, tk)%v(p) = v_new
          species(s)%prtl_tile(ti, tj, tk)%w(p) = w_new

        end do
      end do
      end do
      end do
    end do

    if (mpi_rank .eq. 0) print *, ">>> [8A] imprint_velocity_shear DONE"
  end subroutine imprint_velocity_shear

  subroutine imprint_perturbation()
    implicit none
    integer :: s, ti, tj, tk, p
    integer(kind=2) :: i_idx, k_idx
    real :: x_glob, z_glob, u, v, w
    real :: vx, vy, vz, dvx, dvy
    real :: gamma_old, gamma_new, beta2, scale
    real :: u_new, v_new, w_new
    real :: kz, seed_window
    logical :: perturb_com

    if (.not. enable_perturbation) return

    if (mpi_rank .eq. 0) print *, ">>> [8A] imprint_perturbation..."

    if (abs(B0_ext) > 0.0) then
      kz = 2.0 * PI_F / real(global_mesh%sz)
      perturb_com = .false.
      if (mpi_rank .eq. 0) print *, ">>> [8A] perturbation: localized k=1 sin/cos seed"
    else
      kz = 0.0
      perturb_com = .true.
      if (mpi_rank .eq. 0) print *, ">>> [8A] perturbation: COM dvx kick"
    end if

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp

          i_idx = species(s)%prtl_tile(ti, tj, tk)%xi(p)
          k_idx = species(s)%prtl_tile(ti, tj, tk)%zi(p)
          x_glob = real(i_idx + this_meshblock%ptr%x0) + species(s)%prtl_tile(ti, tj, tk)%dx(p)
          z_glob = real(k_idx + this_meshblock%ptr%z0) + species(s)%prtl_tile(ti, tj, tk)%dz(p)

          u = species(s)%prtl_tile(ti, tj, tk)%u(p)
          v = species(s)%prtl_tile(ti, tj, tk)%v(p)
          w = species(s)%prtl_tile(ti, tj, tk)%w(p)

          gamma_old = sqrt(1.0 + u * u + v * v + w * w)
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

          seed_window = x_seed_window(x_glob)
          dvx = dvx * seed_window
          dvy = dvy * seed_window

          vx = vx + dvx
          vy = vy + dvy

          beta2 = (vx * vx + vy * vy + vz * vz) / (CC * CC)
          if (beta2 >= 1.0) then
            scale = 0.9 / sqrt(beta2)
            vx = vx * scale
            vy = vy * scale
            vz = vz * scale
            beta2 = (vx * vx + vy * vy + vz * vz) / (CC * CC)
          end if

          gamma_new = 1.0 / sqrt(1.0 - beta2)
          u_new = gamma_new * vx / CC
          v_new = gamma_new * vy / CC
          w_new = gamma_new * vz / CC

          species(s)%prtl_tile(ti, tj, tk)%u(p) = u_new
          species(s)%prtl_tile(ti, tj, tk)%v(p) = v_new
          species(s)%prtl_tile(ti, tj, tk)%w(p) = w_new

        end do
      end do
      end do
      end do
    end do

    if (mpi_rank .eq. 0) print *, ">>> [8A] imprint_perturbation DONE"
  end subroutine imprint_perturbation

  subroutine userInitFields()
    implicit none

    ex(:, :, :) = 0.0
    ey(:, :, :) = 0.0
    ez(:, :, :) = 0.0
    bx(:, :, :) = 0.0
    by(:, :, :) = 0.0
    bz(:, :, :) = B0_ext
    jx(:, :, :) = 0.0
    jy(:, :, :) = 0.0
    jz(:, :, :) = 0.0
    if (mpi_rank .eq. 0) then
      print *, ">>> [8A] userInitFields initialized uniform Bz =", B0_ext
    end if
  end subroutine userInitFields

  subroutine userCurrentDeposit(step)
    implicit none
    integer, intent(in), optional :: step
  end subroutine userCurrentDeposit

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
    real :: dt_loc, drive_window
    integer(kind=2) :: temp_i
    real :: temp_r
    real, parameter :: beta2_max = 1.0 - 1.0e-6

    if (.not. use_rotation) return

    dt_loc = 1.0

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp

          i_idx = species(s)%prtl_tile(ti, tj, tk)%xi(p)
          x_glob = real(i_idx + this_meshblock%ptr%x0) + species(s)%prtl_tile(ti, tj, tk)%dx(p)
          x_rel = x_glob - box_xcenter

          u = species(s)%prtl_tile(ti, tj, tk)%u(p)
          v = species(s)%prtl_tile(ti, tj, tk)%v(p)
          w = species(s)%prtl_tile(ti, tj, tk)%w(p)

          gamma_old = sqrt(1.0 + u * u + v * v + w * w)

          vx_old = CC * u / gamma_old
          vy_old = CC * v / gamma_old
          vz_old = CC * w / gamma_old
          vx = vx_old
          vy = vy_old
          vz = vz_old

          drive_window = x_drive_window(x_glob)
          ax = drive_window * (2.0 * Omega * vy + 2.0 * q_shear * Omega * Omega * x_rel)
          ay = drive_window * (-2.0 * Omega * vx)

          vx = vx + ax * dt_loc
          vy = vy + ay * dt_loc

          v2 = vx * vx + vy * vy + vz * vz
          if (v2 >= beta2_max * CC * CC) then
            scale = sqrt((beta2_max * CC * CC) / v2)
            vx = vx * scale
            vy = vy * scale
            vz = vz * scale
            v2 = vx * vx + vy * vy + vz * vz
          end if

          gamma_new = 1.0 / sqrt(1.0 - v2 / (CC * CC))
          species(s)%prtl_tile(ti, tj, tk)%u(p) = gamma_new * vx / CC
          species(s)%prtl_tile(ti, tj, tk)%v(p) = gamma_new * vy / CC
          species(s)%prtl_tile(ti, tj, tk)%w(p) = gamma_new * vz / CC

#if defined(oneD) || defined(twoD) || defined(threeD)
          species(s)%prtl_tile(ti, tj, tk)%dx(p) = species(s)%prtl_tile(ti, tj, tk)%dx(p) + (vx - vx_old) * dt_loc
          temp_i = INT(species(s)%prtl_tile(ti, tj, tk)%dx(p), 2)
          temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti, tj, tk)%dx(p)) + temp_i, REAL(temp_i)) - 1
          temp_i = INT(temp_r, 2)
          species(s)%prtl_tile(ti, tj, tk)%xi(p) = species(s)%prtl_tile(ti, tj, tk)%xi(p) + temp_i
          species(s)%prtl_tile(ti, tj, tk)%dx(p) = species(s)%prtl_tile(ti, tj, tk)%dx(p) - temp_r
#endif

#if defined(twoD) || defined(threeD)
          species(s)%prtl_tile(ti, tj, tk)%dy(p) = species(s)%prtl_tile(ti, tj, tk)%dy(p) + (vy - vy_old) * dt_loc
          temp_i = INT(species(s)%prtl_tile(ti, tj, tk)%dy(p), 2)
          temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti, tj, tk)%dy(p)) + temp_i, REAL(temp_i)) - 1
          temp_i = INT(temp_r, 2)
          species(s)%prtl_tile(ti, tj, tk)%yi(p) = species(s)%prtl_tile(ti, tj, tk)%yi(p) + temp_i
          species(s)%prtl_tile(ti, tj, tk)%dy(p) = species(s)%prtl_tile(ti, tj, tk)%dy(p) - temp_r
#endif

#if defined(threeD)
          species(s)%prtl_tile(ti, tj, tk)%dz(p) = species(s)%prtl_tile(ti, tj, tk)%dz(p) + (vz - vz_old) * dt_loc
          temp_i = INT(species(s)%prtl_tile(ti, tj, tk)%dz(p), 2)
          temp_r = MAX(SIGN(1.0, species(s)%prtl_tile(ti, tj, tk)%dz(p)) + temp_i, REAL(temp_i)) - 1
          temp_i = INT(temp_r, 2)
          species(s)%prtl_tile(ti, tj, tk)%zi(p) = species(s)%prtl_tile(ti, tj, tk)%zi(p) + temp_i
          species(s)%prtl_tile(ti, tj, tk)%dz(p) = species(s)%prtl_tile(ti, tj, tk)%dz(p) - temp_r
#endif

        end do
      end do
      end do
      end do
    end do
  end subroutine userDriveParticles

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, intent(in), optional :: step
    logical, intent(in), optional :: updateE, updateB
    logical :: updateE_, updateB_
    integer :: i, i_glob
    real :: eta, lambda
    real :: x_face, x_edge

    if (present(updateE)) then
      updateE_ = updateE
    else
      updateE_ = .true.
    end if
    if (present(updateB)) then
      updateB_ = updateB
    else
      updateB_ = .true.
    end if

    if ((.not. updateE_) .and. (.not. updateB_)) return

    if ((.not. boundary_reported) .and. (mpi_rank .eq. 0)) then
      print *, ">>> [8A] userFieldBoundaryConditions:"
      print *, "    boundary_x = ", boundary_x, " absorb_x = ", absorb_x, " abs_thick = ", ds_abs
      print *, "    buffer_edge = ", buffer_edge, " kill_edge = ", kill_edge
      boundary_reported = .true.
    end if

    do i = -NGHOST, this_meshblock%ptr%sx - 1 + NGHOST
      i_glob = i + this_meshblock%ptr%x0
      x_edge = real(i_glob)
      x_face = real(i_glob) + 0.5

      if (updateB_) then
        eta = x_buffer_depth(x_edge)
        if (eta > 0.0) then
          lambda = field_damp_kappa * eta**3
          bx(i, :, :) = exp(-lambda) * bx(i, :, :)
        end if

        eta = x_buffer_depth(x_face)
        if (eta > 0.0) then
          lambda = field_damp_kappa * eta**3
          by(i, :, :) = exp(-lambda) * by(i, :, :)
          bz(i, :, :) = exp(-lambda) * bz(i, :, :) + (1.0 - exp(-lambda)) * B0_ext
        end if
      end if

      if (updateE_) then
        eta = x_buffer_depth(x_face)
        if (eta > 0.0) then
          lambda = field_damp_kappa * eta**3
          ex(i, :, :) = exp(-lambda) * ex(i, :, :)
        end if

        eta = x_buffer_depth(x_edge)
        if (eta > 0.0) then
          lambda = field_damp_kappa * eta**3
          ey(i, :, :) = exp(-lambda) * ey(i, :, :)
          ez(i, :, :) = exp(-lambda) * ez(i, :, :)
        end if
      end if
    end do
  end subroutine userFieldBoundaryConditions

  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, intent(in), optional :: step
    integer :: s, ti, tj, tk, p
    real :: x_glob, x_frac, eta, damp

    do s = 1, nspec
      do ti = 1, species(s)%tile_nx
      do tj = 1, species(s)%tile_ny
      do tk = 1, species(s)%tile_nz
        do p = 1, species(s)%prtl_tile(ti, tj, tk)%npart_sp

          if (species(s)%prtl_tile(ti, tj, tk)%proc(p) < 0) cycle

          x_glob = real(species(s)%prtl_tile(ti, tj, tk)%xi(p) + this_meshblock%ptr%x0) + &
                   species(s)%prtl_tile(ti, tj, tk)%dx(p)
          x_frac = x_glob / sx_glob

          if ((x_frac < particle_kill_width_frac) .or. (x_frac > (1.0 - particle_kill_width_frac))) then
            species(s)%prtl_tile(ti, tj, tk)%proc(p) = -1
            cycle
          end if

          eta = x_buffer_depth(x_glob)
          if (eta > 0.0) then
            damp = exp(-particle_damp_kappa * eta**3)
            species(s)%prtl_tile(ti, tj, tk)%u(p) = damp * species(s)%prtl_tile(ti, tj, tk)%u(p)
            species(s)%prtl_tile(ti, tj, tk)%v(p) = damp * species(s)%prtl_tile(ti, tj, tk)%v(p)
            species(s)%prtl_tile(ti, tj, tk)%w(p) = damp * species(s)%prtl_tile(ti, tj, tk)%w(p)
          end if

        end do
      end do
      end do
      end do
    end do
  end subroutine userParticleBoundaryConditions

  subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine readUsrRestart

end module m_userfile
