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

  real :: amplitude
  integer :: mode
  real :: background_T

  private :: amplitude, mode, background_T
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
    return
  end function

  subroutine userReadInput()
    implicit none
    call getInput('problem', 'amplitude', amplitude)
    call getInput('problem', 'mode', mode, 1)
    call getInput('problem', 'background_T', background_T, 0.0)
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3

    return
  end function

  subroutine userInitParticles()
    implicit none
    integer :: s, ti, tj, tk, p
    real :: kx, phase, vA, delta_beta
    real :: xp
    type(region) :: back_region
    procedure(spatialDistribution), pointer :: spat_distr_ptr => null()

    spat_distr_ptr => userSpatialDistribution

    back_region % x_min = 0.0
    back_region % x_max = global_mesh % sx
#ifdef twoD
    back_region % y_min = 0.0
    back_region % y_max = global_mesh % sy
#endif

    ! uniform pair plasma with no initial current
    call fillRegionWithThermalPlasma(back_region, (/1, 2/), 2, ppc0, &
                                     background_T, zero_current=.true., &
                                     spat_distr_ptr=spat_distr_ptr)

    kx = 2.0 * M_PI * REAL(mode) / REAL(global_mesh % sx)
    vA = CC * sqrt(sigma / (1.0 + sigma))
    delta_beta = amplitude * vA / CC

    ! small y-velocity perturbation to seed the shear Alfvén mode
    do s = 1, nspec
      do ti = 1, species(s) % tile_nx
        do tj = 1, species(s) % tile_ny
          do tk = 1, species(s) % tile_nz
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              xp = REAL(species(s) % prtl_tile(ti, tj, tk) % xi(p)) + &
                   species(s) % prtl_tile(ti, tj, tk) % dx(p)
              xp = xp + REAL(this_meshblock % ptr % x0)
              phase = kx * xp
              species(s) % prtl_tile(ti, tj, tk) % v(p) = &
                   species(s) % prtl_tile(ti, tj, tk) % v(p) + delta_beta * sin(phase)
            end do
          end do
        end do
      end do
    end do
  end subroutine userInitParticles

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    integer :: i_glob, j_glob, k_glob
    real :: kx, phase, vA, deltaB

    ex(:, :, :) = 0.0; ey(:, :, :) = 0.0; ez(:, :, :) = 0.0
    bx(:, :, :) = 1.0; by(:, :, :) = 0.0; bz(:, :, :) = 0.0
    jx(:, :, :) = 0.0; jy(:, :, :) = 0.0; jz(:, :, :) = 0.0

    kx = 2.0 * M_PI * REAL(mode) / REAL(global_mesh % sx)
    vA = CC * sqrt(sigma / (1.0 + sigma))

    do i = 0, this_meshblock % ptr % sx - 1
      i_glob = i + this_meshblock % ptr % x0
      phase = kx * REAL(i_glob)
      deltaB = amplitude * sin(phase)
      do j = 0, this_meshblock % ptr % sy - 1
        j_glob = j + this_meshblock % ptr % y0
        do k = 0, this_meshblock % ptr % sz - 1
          k_glob = k + this_meshblock % ptr % z0
          by(i, j, k) = deltaB
          ez(i, j, k) = - (vA / CC) * amplitude * &
                        sin(kx * (REAL(i_glob) - 0.5))
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
