module m_currentdeposit
  use m_globalnamespace
  use m_aux
  use m_domain
  use m_fields
  use m_particles
  use m_errors, only: throwError
  implicit none

  logical, save :: deposit_trace_enable = .false.
  logical, save :: deposit_trace_init = .false.
contains
  subroutine maybeEnableDepositTrace()
    implicit none
    character(len=32) :: env
    integer :: stat
    if (deposit_trace_init) return
    deposit_trace_init = .true.
    env = ''
    call get_environment_variable('TRISTAN_DEPOSIT_TRACE', env, status=stat)
    if (stat == 0) then
      if (len_trim(env) > 0 .and. trim(env) /= '0') deposit_trace_enable = .true.
    end if
  end subroutine maybeEnableDepositTrace

  subroutine resetCurrents()
    implicit none
    jx(:, :, :) = 0; jy(:, :, :) = 0; jz(:, :, :) = 0
  end subroutine resetCurrents

  subroutine depositCurrents()
    ! DEP_PRT [particle-dependent]
    implicit none
    integer :: s, p, ti, tj, tk
    integer(kind=2), pointer, contiguous :: pt_xi(:), pt_yi(:), pt_zi(:)
    real, pointer, contiguous :: pt_dx(:), pt_dy(:), pt_dz(:), pt_wei(:)
    real :: xr, yr, zr, x1, y1, z1, x2, y2, z2
    real :: gamma_inv, temp_charge, weighted_charge
    integer(kind=2) :: i1, i2, j1, j2, k1, k2
    integer(kind=2) :: i1p1, i2p1, j1p1, j2p1, k1p1, k2p1
    real :: Wx1, Wy1, Wz1, Wx2, Wy2, Wz2
    real :: onemWx1, onemWy1, onemWz1, onemWx2, onemWy2, onemWz2
    real :: Fx1, Fy1, Fz1, Fx2, Fy2, Fz2
    real, pointer, contiguous :: pt_u(:), pt_v(:), pt_w(:)
    integer :: ilb, iub, jlb, jub, klb, kub

    call maybeEnableDepositTrace()

    ilb = lbound(jx, 1); iub = ubound(jx, 1)
    jlb = lbound(jx, 2); jub = ubound(jx, 2)
    klb = lbound(jx, 3); kub = ubound(jx, 3)

    do s = 1, nspec ! loop over species
      if ((species(s) % ch_sp .eq. 0) .or. (.not. species(s) % deposit_sp)) cycle
      do tk = 1, species(s) % tile_nz
        do tj = 1, species(s) % tile_ny
          do ti = 1, species(s) % tile_nx
            pt_xi => species(s) % prtl_tile(ti, tj, tk) % xi
            pt_yi => species(s) % prtl_tile(ti, tj, tk) % yi
            pt_zi => species(s) % prtl_tile(ti, tj, tk) % zi

            pt_dx => species(s) % prtl_tile(ti, tj, tk) % dx
            pt_dy => species(s) % prtl_tile(ti, tj, tk) % dy
            pt_dz => species(s) % prtl_tile(ti, tj, tk) % dz

            pt_u => species(s) % prtl_tile(ti, tj, tk) % u
            pt_v => species(s) % prtl_tile(ti, tj, tk) % v
            pt_w => species(s) % prtl_tile(ti, tj, tk) % w
            pt_wei => species(s) % prtl_tile(ti, tj, tk) % weight

            temp_charge = species(s) % ch_sp * unit_ch / B_norm
            do p = 1, species(s) % prtl_tile(ti, tj, tk) % npart_sp
              ! push the particle back
              gamma_inv = 1.0 / sqrt(1.0 + pt_u(p)**2 + pt_v(p)**2 + pt_w(p)**2)

              x2 = REAL(pt_xi(p)) + pt_dx(p); y2 = REAL(pt_yi(p)) + pt_dy(p); z2 = REAL(pt_zi(p)) + pt_dz(p)
              x1 = x2 - pt_u(p) * CC * gamma_inv; y1 = y2 - pt_v(p) * CC * gamma_inv; z1 = z2 - pt_w(p) * CC * gamma_inv

#ifdef oneD
              i1 = INT(FLOOR(x1), 2); i2 = pt_xi(p)
              j1 = 0; j2 = 0
              k1 = 0; k2 = 0
              i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
#elif defined(twoD)
              i1 = INT(FLOOR(x1), 2); i2 = pt_xi(p)
              j1 = INT(FLOOR(y1), 2); j2 = pt_yi(p)
              k1 = 0; k2 = 0
              i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
              j1p1 = j1 + 1_2; j2p1 = j2 + 1_2
#elif defined(threeD)
              i1 = INT(FLOOR(x1), 2); i2 = pt_xi(p)
              j1 = INT(FLOOR(y1), 2); j2 = pt_yi(p)
              k1 = INT(FLOOR(z1), 2); k2 = pt_zi(p)
              i1p1 = i1 + 1_2; i2p1 = i2 + 1_2
              j1p1 = j1 + 1_2; j2p1 = j2 + 1_2
              k1p1 = k1 + 1_2; k2p1 = k2 + 1_2
#endif

              weighted_charge = pt_wei(p) * temp_charge

              if (deposit_trace_enable) then
                if ((x1 /= x1) .or. (x2 /= x2) .or. (y1 /= y1) .or. (y2 /= y2) .or. (z1 /= z1) .or. (z2 /= z2)) then
                  print *, "DEPOSIT_TRACE: NaN in trajectory coords", " s=", s, " p=", p, &
                           " x1/x2=", x1, x2, " y1/y2=", y1, y2, " z1/z2=", z1, z2
                  call throwError("DEPOSIT_TRACE: NaN in x1/x2/y1/y2/z1/z2")
                end if

                if ((i1 < ilb) .or. (i1 > iub) .or. (i1p1 < ilb) .or. (i1p1 > iub) .or. &
                    (i2 < ilb) .or. (i2 > iub) .or. (i2p1 < ilb) .or. (i2p1 > iub) .or. &
                    (j1 < jlb) .or. (j1 > jub) .or. (j1p1 < jlb) .or. (j1p1 > jub) .or. &
                    (j2 < jlb) .or. (j2 > jub) .or. (j2p1 < jlb) .or. (j2p1 > jub) .or. &
                    (k1 < klb) .or. (k1 > kub) .or. (k2 < klb) .or. (k2 > kub)) then
                  print *, "DEPOSIT_TRACE: OOB indices", " s=", s, " tile=", ti, tj, tk, " p=", p, &
                           " i1/i2=", i1, i2, " j1/j2=", j1, j2, " k1/k2=", k1, k2, &
                           " bounds i/j/k=", ilb, iub, jlb, jub, klb, kub, &
                           " x1/x2=", x1, x2, " y1/y2=", y1, y2, " z1/z2=", z1, z2
                  call throwError("DEPOSIT_TRACE: out-of-bounds deposit indices")
                end if
              end if

! this "function" takes
! ... the start and end coordinates: `x1`, `x2`, `y1`, `y2`, `z1`, `z2` ...
! ... the start and end cells: `i1`, `i2`, `j1`, `j2`, `k1`, `k2` ...
! ... the start and end cells + 1: `i1p1`, `i2p1` etc ...
! ... the weighted_chargeed charge: `weighted_charge = weight * charge_sp * unit_charge / Bnorm`
! ... and deposits proper currents to corresponding components
#include "zigzag_deposit.F08"
            end do
            pt_xi => null(); pt_yi => null(); pt_zi => null()
            pt_dx => null(); pt_dy => null(); pt_dz => null()
            pt_wei => null()
            pt_u => null(); pt_v => null(); pt_w => null()
          end do
        end do
      end do
    end do ! species loop

    call printDiag("depositCurrents()", 2)
  end subroutine depositCurrents
end module m_currentdeposit
