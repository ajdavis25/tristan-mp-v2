module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_thermalplasma
  use m_particlelogistics
  implicit none

  real :: B0_harris, a_harris, n0_sheet, n_bg_harris
  real :: Ti_harris, Te_harris
  real :: sheet_center

  private :: B0_harris, a_harris, n0_sheet, n_bg_harris
  private :: Ti_harris, Te_harris, sheet_center
  private :: userSpatialDistribution

contains

  subroutine userDeallocate()
    implicit none
  end subroutine userDeallocate

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    userSLBload = 0.0
    return
  end function userSLBload

  subroutine userReadInput()
    implicit none
    call getInput('problem', 'B0',  B0_harris)
    call getInput('problem', 'a',   a_harris)
    call getInput('problem', 'n0',  n0_sheet)
    call getInput('problem', 'n_bg', n_bg_harris, 0.0)
    call getInput('problem', 'Ti',  Ti_harris, 0.0)
    call getInput('problem', 'Te',  Te_harris, 0.0)

    if (a_harris .le. 0.0) then
      call throwError('ERROR: Harris half-thickness `a` must be > 0.')
    end if
    if (n0_sheet .lt. 0.0) then
      call throwError('ERROR: Harris peak density `n0` must be >= 0.')
    end if
    if (n_bg_harris .lt. 0.0) then
      call throwError('ERROR: Harris background density `n_bg` must be >= 0.')
    end if
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    userSpatialDistribution = 1.0
    if (present(y_glob)) then
      ! Harris density ∝ sech^2((y - y0)/a)
      userSpatialDistribution = 1.0 / cosh((y_glob - sheet_center) / a_harris)**2
    end if
  end function userSpatialDistribution

  subroutine userInitParticles()
    implicit none
    real :: base_ppc, bg_density, sheet_density
    real :: shift_beta_raw, shift_beta, shift_gamma
    integer :: drift_dir
    type(region) :: full_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    sheet_center = 0.5 * REAL(global_mesh % sy)
    spat_distr_ptr => userSpatialDistribution

    full_region % x_min = 0.0
    full_region % y_min = 0.0
    full_region % x_max = REAL(global_mesh % sx)
    full_region % y_max = REAL(global_mesh % sy)
#ifdef threeD
    full_region % z_min = 0.0
    full_region % z_max = REAL(global_mesh % sz)
#endif

    ! base_ppc is per species in this setup
    base_ppc = 0.5 * ppc0
    bg_density    = base_ppc * n_bg_harris
    sheet_density = base_ppc * n0_sheet

    ! Background pair plasma (uniform, zero net current)
    if (bg_density .gt. 0.0) then
      call fillRegionWithThermalPlasma(full_region, (/1/), 1, bg_density, Te_harris, zero_current=.true.)
      call fillRegionWithThermalPlasma(full_region, (/2/), 1, bg_density, Ti_harris, zero_current=.true.)
    end if

    ! Harris sheet drifting population
    if (sheet_density .gt. 0.0) then
      ! Drift beta chosen to balance B0, n0, and sheet thickness a; clamp to ~0.75–0.9 to drive tearing faster.
      shift_beta_raw = (abs(B0_harris) * sqrt(sigma) * c_omp) / (a_harris * n0_sheet)
      shift_beta = max(0.75, min(0.9, shift_beta_raw))
      if (shift_beta_raw .ge. 1.0) then
        call throwError('ERROR: Harris drift `shift_beta` >= 1. Adjust B0, a, or n0.')
      end if
      shift_gamma = 1.0 / sqrt(1.0 - shift_beta**2)

      ! Drift along ±z depending on sign of B0
      drift_dir = -3 * INT(sign(1.0, B0_harris))

      call fillRegionWithThermalPlasma(full_region, (/1/), 1, sheet_density, Te_harris, &
                                       shift_gamma=shift_gamma, shift_dir=drift_dir, &
                                       spat_distr_ptr=spat_distr_ptr)
      call fillRegionWithThermalPlasma(full_region, (/2/), 1, sheet_density, Ti_harris, &
                                       shift_gamma=shift_gamma, shift_dir=drift_dir, &
                                       spat_distr_ptr=spat_distr_ptr)
    end if
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k, kmin, kmax
    real :: x_glob, y_glob, sy_glob, profile
    real :: delta, kx, ky, az
    real, parameter :: pi = 3.141592653589793

    ex(:, :, :) = 0.0; ey(:, :, :) = 0.0; ez(:, :, :) = 0.0
    bx(:, :, :) = 0.0; by(:, :, :) = 0.0; bz(:, :, :) = 0.0
    jx(:, :, :) = 0.0; jy(:, :, :) = 0.0; jz(:, :, :) = 0.0

    sy_glob = REAL(global_mesh % sy)
    sheet_center = 0.5 * sy_glob

    ! Harris equilibrium Bx(y) = B0 tanh((y - y0)/a)
    kmin = 0; kmax = 0
#ifdef threeD
    kmin = -NGHOST
    kmax = this_meshblock % ptr % sz - 1 + NGHOST
#endif

    do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST
      y_glob = REAL(j + this_meshblock % ptr % y0) + 0.5
      profile = tanh((y_glob - sheet_center) / a_harris)
      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
        do k = kmin, kmax
          bx(i, j, k) = B0_harris * profile
        end do
      end do
    end do

    ! Seed a long-wavelength tearing perturbation via δA_z
    delta = 0.10 * B0_harris
    kx = 2.0 * pi / REAL(global_mesh % sx)
    ky = pi / REAL(global_mesh % sy)

    do j = -NGHOST, this_meshblock % ptr % sy - 1 + NGHOST
      y_glob = REAL(j + this_meshblock % ptr % y0) + 0.5
      do i = -NGHOST, this_meshblock % ptr % sx - 1 + NGHOST
        x_glob = REAL(i + this_meshblock % ptr % x0) + 0.5
        do k = kmin, kmax
          az = delta * cos(kx * x_glob) * cos(ky * (y_glob - sheet_center))
          bx(i, j, k) = bx(i, j, k) + az * ky
          by(i, j, k) = by(i, j, k) - az * kx
        end do
      end do
    end do
  end subroutine userInitFields

  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = 0.0; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields

  subroutine writeUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine writeUsrRestart

  subroutine readUsrRestart(rst_file)
    implicit none
    integer, intent(in) :: rst_file
  end subroutine readUsrRestart

end module m_userfile
