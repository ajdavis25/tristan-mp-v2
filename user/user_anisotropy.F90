module m_userfile
  use m_globalnamespace
  use m_aux
  use m_readinput
  use m_domain
  use m_particles
  use m_fields
  use m_helpers
  use m_particlelogistics
  use m_thermalplasma
  use m_errors
  implicit none

  !--- PRIVATE variables -----------------------------------------!
  real :: tpar0, tperp0, b0_aniso, n0_uniform, noise_frac, noise_phase
  private :: tpar0, tperp0, b0_aniso, n0_uniform, noise_frac, noise_phase
  !...............................................................!

  !--- PRIVATE functions -----------------------------------------!
  private :: userSpatialDistribution
  private :: loadBiMaxwellian
  !...............................................................!
contains
  !--- initialization -----------------------------------------!
  subroutine userReadInput()
    implicit none
    real :: b2_phys, firehose_rhs, beta_par, mirror_rhs
    character(len=128) :: msg
    call getInput('problem', 'Tpar0', tpar0)
    call getInput('problem', 'Tperp0', tperp0)
    call getInput('problem', 'B0', b0_aniso)
    call getInput('problem', 'n0', n0_uniform, 1.0)
    call getInput('problem', 'noise_frac', noise_frac, 0.0)
    noise_phase = 2.0 * M_PI * random(dseed)

    b2_phys = (b0_aniso * B_norm)**2
    firehose_rhs = 0.0
    beta_par = HUGE(1.0)
    mirror_rhs = 0.0
    if (b2_phys .gt. 0.0) then
      firehose_rhs = b2_phys / (4.0 * M_PI * n0_uniform)
      beta_par = 8.0 * M_PI * n0_uniform * tpar0 / b2_phys
      mirror_rhs = 1.0 / beta_par
    end if
    if (mpi_rank .eq. 0) then
      write (msg, '(A,F10.4,A,F10.4)') 'Firehose LHS (Tpar-Tperp)=', tpar0 - tperp0, &
        '  RHS=B^2/(4πn)=', firehose_rhs
      call printDiag(trim(msg), 1)
      write (msg, '(A,F10.4,A,F10.4,A,F10.4)') 'Mirror LHS (Tperp/Tpar-1)=', &
        (tperp0 / tpar0 - 1.0), '  beta_par=', beta_par, '  RHS=1/beta_par=', mirror_rhs
      call printDiag(trim(msg), 1)
    end if
  end subroutine userReadInput

  function userSpatialDistribution(x_glob, y_glob, z_glob, &
                                   dummy1, dummy2, dummy3)
    real :: userSpatialDistribution
    real, intent(in), optional :: x_glob, y_glob, z_glob
    real, intent(in), optional :: dummy1, dummy2, dummy3
    userSpatialDistribution = 1.0
    return
  end function

  function userSLBload(x_glob, y_glob, z_glob, &
                       dummy1, dummy2, dummy3)
    real :: userSLBload
    ! global coordinates
    real, intent(in), optional :: x_glob, y_glob, z_glob
    ! global box dimensions
    real, intent(in), optional :: dummy1, dummy2, dummy3
    userSLBload = 1.0
    return
  end function

  subroutine userInitParticles()
    implicit none
    type(region) :: full_region
    real :: ndens_per_species
    integer :: s
    integer, allocatable :: fill_species(:)

    full_region % x_min = 0.0
    full_region % x_max = REAL(global_mesh % sx)
#ifdef twoD
    full_region % y_min = 0.0
    full_region % y_max = REAL(global_mesh % sy)
#endif
#ifdef threeD
    full_region % z_min = 0.0
    full_region % z_max = REAL(global_mesh % sz)
#endif

    ndens_per_species = n0_uniform * ppc0
    allocate (fill_species(nspec))
    do s = 1, nspec
      fill_species(s) = s
    end do
    call loadBiMaxwellian(full_region, fill_species, nspec, ndens_per_species, &
                          tpar0, tperp0)
    deallocate (fill_species)
  end subroutine userInitParticles

  subroutine loadBiMaxwellian(fill_region, fill_species, num_species, ndens_sp, tpar, tperp)
    implicit none
    type(region), intent(in) :: fill_region
    integer, intent(in) :: num_species
    integer, intent(in) :: fill_species(num_species)
    real, intent(in) :: ndens_sp, tpar, tperp

    real :: fill_xmin, fill_xmax, &
            fill_ymin, fill_ymax, &
            fill_zmin, fill_zmax
    real :: num_part_r
    integer :: num_part, n, s, spec_
    integer(kind=2) :: xi_, yi_, zi_
    real :: u_par, v_perp, w_perp
    real :: x_, y_, z_, dx_, dy_, dz_, dummy1, dummy2
    type(maxwellian) :: maxw_par, maxw_perp

    maxw_par % dimension = 1
    maxw_par % shift_flag = .false.
    maxw_par % generated = .false.

    maxw_perp % dimension = 1
    maxw_perp % shift_flag = .false.
    maxw_perp % generated = .false.

#ifdef oneD
    call globalToLocalCoords(fill_region % x_min, 0.0, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, 0.0, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin)
#elif defined(twoD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, 0.0, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, 0.0, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin) * (fill_ymax - fill_ymin)
#elif defined(threeD)
    call globalToLocalCoords(fill_region % x_min, fill_region % y_min, fill_region % z_min, &
                             fill_xmin, fill_ymin, fill_zmin, adjustQ=.true.)
    call globalToLocalCoords(fill_region % x_max, fill_region % y_max, fill_region % z_max, &
                             fill_xmax, fill_ymax, fill_zmax, adjustQ=.true.)
    num_part_r = REAL(ndens_sp) * (fill_xmax - fill_xmin) * &
                 (fill_ymax - fill_ymin) * (fill_zmax - fill_zmin)
#endif

    if (num_part_r .lt. 10.0) then
      if (num_part_r .ne. 0.0) then
        num_part_r = poisson(num_part_r)
      else
        num_part_r = 0.0
      end if
    else
      num_part_r = CEILING(num_part_r)
    end if
    num_part = INT(num_part_r)

    n = 0
    do while (n .lt. num_part)
      call generateCoordInRegion(fill_xmin, fill_xmax, fill_ymin, fill_ymax, fill_zmin, fill_zmax, &
                                 x_, y_, z_, xi_, yi_, zi_, dx_, dy_, dz_)
      do s = 1, num_species
        spec_ = fill_species(s)
        if ((spec_ .le. 0) .or. (spec_ .gt. nspec)) then
          call throwError('Wrong species specified in loadBiMaxwellian.')
        end if
        maxw_par % temperature = tpar / species(spec_) % m_sp
        maxw_perp % temperature = tperp / species(spec_) % m_sp
        call generateFromMaxwellian(maxw_par, u_par, dummy1, dummy2)
        call generateFromMaxwellian(maxw_perp, v_perp, dummy1, dummy2)
        call generateFromMaxwellian(maxw_perp, w_perp, dummy1, dummy2)
        call createParticle(spec_, xi_, yi_, zi_, dx_, dy_, dz_, u_par, v_perp, w_perp)
      end do
      n = n + 1
    end do
  end subroutine loadBiMaxwellian

  subroutine userInitFields()
    implicit none
    integer :: i, j, k
    real :: y_glob, ky, noise_amp, phase

    ex(:, :, :) = 0.0; ey(:, :, :) = 0.0; ez(:, :, :) = 0.0
    bx(:, :, :) = 0.0; by(:, :, :) = 0.0; bz(:, :, :) = 0.0
    jx(:, :, :) = 0.0; jy(:, :, :) = 0.0; jz(:, :, :) = 0.0

    ky = 0.0
    if (global_mesh % sy .gt. 0) ky = 2.0 * M_PI / MAX(1.0, REAL(global_mesh % sy))
    noise_amp = noise_frac * b0_aniso

    do j = this_meshblock % ptr % j1, this_meshblock % ptr % j2
      y_glob = REAL(j + this_meshblock % ptr % y0) + 0.5
      do i = this_meshblock % ptr % i1, this_meshblock % ptr % i2
        do k = this_meshblock % ptr % k1, this_meshblock % ptr % k2
          bx(i, j, k) = b0_aniso
          by(i, j, k) = 0.0
          bz(i, j, k) = 0.0
          if (noise_amp .ne. 0.0) then
            phase = ky * y_glob + noise_phase
            bz(i, j, k) = noise_amp * sin(phase)
          end if
        end do
      end do
    end do
  end subroutine userInitFields
  !............................................................!

  !--- driving ------------------------------------------------!
  subroutine userCurrentDeposit(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userCurrentDeposit

  subroutine userDriveParticles(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userDriveParticles

  subroutine userExternalFields(xp, yp, zp, &
                                ex_ext, ey_ext, ez_ext, &
                                bx_ext, by_ext, bz_ext)
    implicit none
    real, intent(in) :: xp, yp, zp
    real, intent(out) :: ex_ext, ey_ext, ez_ext
    real, intent(out) :: bx_ext, by_ext, bz_ext
    ex_ext = 0.0; ey_ext = 0.0; ez_ext = 0.0
    bx_ext = b0_aniso; by_ext = 0.0; bz_ext = 0.0
  end subroutine userExternalFields
  !............................................................!

  !--- boundaries ---------------------------------------------!
  subroutine userParticleBoundaryConditions(step)
    implicit none
    integer, optional, intent(in) :: step
  end subroutine userParticleBoundaryConditions

  subroutine userFieldBoundaryConditions(step, updateE, updateB)
    implicit none
    integer, optional, intent(in) :: step
    logical, optional, intent(in) :: updateE, updateB
  end subroutine userFieldBoundaryConditions
  !............................................................!

#include "optional.F"
end module m_userfile
