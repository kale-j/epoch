! Copyright (C) 2009-2019 University of Warwick
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

MODULE analytic_pulse

#ifdef APT_VACUUM
  USE shared_data

  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_analytic_pulse(boundary, analytic_pulse)

    INTEGER, INTENT(IN) :: boundary
    TYPE(analytic_pulse_block), INTENT(INOUT) :: analytic_pulse

    analytic_pulse%boundary = boundary
    analytic_pulse%omega = -1.0_num
    analytic_pulse%a0 = -1.0_num
    analytic_pulse%w0 = -1.0_num
    analytic_pulse%tau = -1.0_num
    analytic_pulse%bf = 0.0_num
    analytic_pulse%xf = 0.0_num
    analytic_pulse%yf = 0.0_num
    analytic_pulse%zf = 0.0_num
    analytic_pulse%t0 = 0.0_num
    analytic_pulse%ph0 = 0.0_num
    analytic_pulse%sg = 1
    NULLIFY(analytic_pulse%next)

  END SUBROUTINE init_analytic_pulse


  SUBROUTINE deallocate_analytic_pulse(analytic_pulse)

    TYPE(analytic_pulse_block), POINTER :: analytic_pulse
    DEALLOCATE(analytic_pulse)

  END SUBROUTINE deallocate_analytic_pulse



  SUBROUTINE deallocate_analytic_pulses

    TYPE(analytic_pulse_block), POINTER :: current, next

    current => analytic_pulses
    DO WHILE(ASSOCIATED(current))
      next => current%next
      CALL deallocate_analytic_pulse(current)
      current => next
    END DO

  END SUBROUTINE deallocate_analytic_pulses


  ! Subroutine to attach a created analytic_pulse object to the correct boundary
  SUBROUTINE attach_analytic_pulse(analytic_pulse)

    TYPE(analytic_pulse_block), POINTER :: analytic_pulse
    TYPE(analytic_pulse_block), POINTER :: current
    INTEGER :: boundary

    boundary = analytic_pulse%boundary

    n_analytic_pulses(boundary) = n_analytic_pulses(boundary) + 1

    IF (ASSOCIATED(analytic_pulses)) THEN
      current => analytic_pulses
      DO WHILE(ASSOCIATED(current%next))
        current => current%next
      END DO
      current%next => analytic_pulse
    ELSE
      analytic_pulses => analytic_pulse
    END IF

  END SUBROUTINE attach_analytic_pulse


  SUBROUTINE setup_analytic_pulses
    ! called after dt is set
    REAL(num) :: f, rfac
    TYPE(analytic_pulse_block), POINTER :: current
    f = c*dt/dx

    ! only x_min for now
    current => analytic_pulses
    DO WHILE(ASSOCIATED(current))
      IF (current%boundary == c_bd_x_min) THEN
        rfac = 0.25*(current%omega*dx/c)**2
        ! b0 set from a0
        current%b0 = current%a0 * current%omega * m0/q0
        ! e0 with numerical correction
        current%e0 = current%b0 * c * (1 + 0.5*f**2*rfac)
        ! e1, b1 from Lortentz transform
        current%e1 = current%e0 * (1-current%bf)
        current%b1 = current%b0 * (1-current%bf)
        current%rayl = 0.5*(1-current%bf)*current%w0**2*current%omega/c
        ! numerical phase and group velocity
        current%vph = c*(1-(1-f**2)*rfac/6)
        current%vg = c*(1-0.5*(1-f**2)*rfac)
        ! xf adjusted so focus will be at input xf at time t0
        current%xf = (1-current%bf)*current%xf - c*current%bf*current%t0 &
             + 0.5*c*current%bf*dt + current%bf*x_min
        ! t0 correction for evaluation time: t0 -> t0 - 0.5*dt
        ! adjusted so that intensity peak is at input xf at time t0
        current%t0 = current%t0 - 0.5*dt - x_min/current%vg
        current%ph0 = current%ph0 + 0.5*current%omega*dt + current%omega*x_min/current%vph
      END IF
      current => current%next
    END DO

  END SUBROUTINE setup_analytic_pulses


  SUBROUTINE analytic_pulse_total_fields

    CALL analytic_pulse_update_arrays

    ex_total = ex_total + ex
    ey_total = ey_total + ey
    ez_total = ez_total + ez
    bx_total = bx_total + bx
    by_total = by_total + by
    bz_total = bz_total + bz

  END SUBROUTINE analytic_pulse_total_fields


  SUBROUTINE analytic_pulse_update_arrays

    TYPE(analytic_pulse_block), POINTER :: current
    TYPE(parameter_pack) :: parameters
    INTEGER :: n,i,j,k,err
    ! copy for convenience
    REAL(num) :: iw02, rayl, irayl, irayl2, tenv_norm
    REAL(num) :: itau_ivg, xf, yf, zf, omega_ivph, pht
    REAL(num) :: e0_amp, b0_amp, e1_amp, b1_amp
    INTEGER :: sg
    REAL(num) :: amp, phase

    ex_total = 0
    ey_total = 0
    ez_total = 0
    bx_total = 0
    by_total = 0
    bz_total = 0

    current => analytic_pulses

    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the analytic pulse
      ! loop over position
      SELECT CASE(current%boundary)

        CASE(c_bd_x_min)
          iw02 = 1.0_num/(current%w0*current%w0)
          rayl = current%rayl
          irayl = 1.0_num/current%rayl
          irayl2 = 1.0_num/(current%rayl*current%rayl)
          tenv_norm = (time-current%t0)/current%tau
          itau_ivg = 1.0_num/(current%tau*current%vg)
          sg = 2*current%sg
          xf = current%xf + c*current%bf*time
          yf = current%yf
          zf = current%zf
          omega_ivph = current%omega/current%vph
          pht = current%ph0 + current%omega*time
          e0_amp = current%e0
          b0_amp = current%b0
          e1_amp = current%e1
          b1_amp = current%b1

          ! ex (two parts)
          DO i = 1-ng,nx+ng
             amp = -e1_amp * irayl * yf &
                  * EXP( -iw02*(yf*yf+zf*zf)/(1+irayl2*(xf-x(i))**2) &
                  -(tenv_norm-itau_ivg*x(i))**sg ) / ((1+irayl2*(xf-x(i))**2)*SQRT(1+irayl2*(xf-x(i))**2))
             phase = pht-omega_ivph*x(i) - ATAN2(xf-x(i),rayl) &
                  + iw02*irayl*(yf*yf+zf*zf)*(xf-x(i))/(1+irayl2*(xf-x(i))**2)
             ex_total(i) = ex_total(i) + amp*( irayl*(xf-x(i))*SIN(phase) + COS(phase) )
          END DO

          ! ey
          DO i = 1-ng,nx+ng
             ey_total(i) = ey_total(i) + e0_amp * &
                  EXP( -iw02*(yf*yf+zf*zf)/(1+irayl2*(xf-xb(i))**2) &
                  -(tenv_norm-itau_ivg*xb(i))**sg ) / SQRT(1+irayl2*(xf-xb(i))**2) &
                  *SIN( pht-omega_ivph*xb(i) - ATAN2(xf-xb(i),rayl) &
                  + iw02*irayl*(yf*yf+zf*zf)*(xf-xb(i))/(1+irayl2*(xf-xb(i))**2) )
          END DO

          ! bx (two parts)
          DO i = 1-ng,nx+ng
             amp = -b1_amp * irayl * zf &
                  * EXP( -iw02*(yf*yf+zf*zf)/(1+irayl2*(xf-xb(i))**2) &
                  -(tenv_norm-itau_ivg*xb(i))**sg ) / ((1+irayl2*(xf-xb(i))**2)*SQRT(1+irayl2*(xf-xb(i))**2))
             phase = pht-omega_ivph*xb(i) - ATAN2(xf-xb(i),rayl) &
                  + iw02*irayl*(yf*yf+zf*zf)*(xf-xb(i))/(1+irayl2*(xf-xb(i))**2)
             bx_total(i) = bx_total(i) + amp*( irayl*(xf-xb(i))*SIN(phase) + COS(phase) )
          END DO

          ! bz
          DO i = 1-ng,nx+ng
             bz_total(i) = bz_total(i) + b0_amp * &
                  EXP( -iw02*(yf*yf+zf*zf)/(1+irayl2*(xf-x(i))**2) &
                  -(tenv_norm-itau_ivg*x(i))**sg ) / SQRT(1+irayl2*(xf-x(i))**2) &
                  *SIN( pht-omega_ivph*x(i) - ATAN2(xf-x(i),rayl) &
                  + iw02*irayl*(yf*yf+zf*zf)*(xf-x(i))/(1+irayl2*(xf-x(i))**2) )
          END DO


        END SELECT
        current => current%next
      END DO

  END SUBROUTINE analytic_pulse_update_arrays

#endif

END MODULE analytic_pulse
