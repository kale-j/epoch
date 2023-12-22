! Copyright (C) 2023 Kale Weichman
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

#if defined(APT_VACUUM) || defined(APT_PLASMA)
  USE shared_data
  
  IMPLICIT NONE

CONTAINS

  SUBROUTINE init_analytic_pulse(boundary, analytic_pulse)

    INTEGER, INTENT(IN) :: boundary
    TYPE(analytic_pulse_block), INTENT(INOUT) :: analytic_pulse

    analytic_pulse%boundary = boundary
    analytic_pulse%omega = -1.0_num
    analytic_pulse%a0 = -1.0_num
    analytic_pulse%tau = -1.0_num
    analytic_pulse%t0 = 0.0_num
    analytic_pulse%ph0 = 0.0_num
    analytic_pulse%sg = 1
#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)
    analytic_pulse%w0 = -1.0_num
    analytic_pulse%bf = 0.0_num
    analytic_pulse%xf = 0.0_num
    analytic_pulse%yf = 0.0_num
#ifdef APT_VACUUM_GAUSS    
    analytic_pulse%zf = 0.0_num
#endif
#endif
#ifdef APT_PLASMA
    analytic_pulse%wp2norm = -1.0_num
#endif
    
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
    REAL(num) :: f2, rfac, wp2n, cvph2, s0
    TYPE(analytic_pulse_block), POINTER :: current
    
    ! only x_min for now
    current => analytic_pulses
    DO WHILE(ASSOCIATED(current))
      IF (current%boundary == c_bd_x_min) THEN
        f2 = (c*dt/dx)**2
        rfac = 0.25*(current%omega*dx/c)**2

        ! e0 set from a0
        current%e0 = current%a0 * current%omega * m0 * c/q0
        ! phase velocity, group velocity, and b0, from dispersion relation
#ifdef APT_VACUUM
        wp2n = 0
#else
        ! plasma
        wp2n = current%wp2norm
#endif
        ! general formulation for phase velocity based on dispersion relation has a parameter related to particle shape
        ! s0 = 0.5*(b+f2*(a-1)); for triangle particle shape, b=2 and a=1
        ! using the default value with a different particle shape may impair the stability of the simulation in some cases
        s0 = 1.0_num
        current%vph = c*(1-(1-f2 + wp2n*s0/(1-wp2n))*rfac/6)/SQRT(1-wp2n)
        current%vg = c*SQRT(1-wp2n)*(1-(3*(1-f2)-wp2n*(2-2*f2-s0)-wp2n**2*s0/(1-wp2n))*rfac/6)
        cvph2 = (c/current%vph)**2
        ! b0 from e0 and phase velocity
        current%b0 = current%e0 * (1-(cvph2+2*f2)*rfac/6)/current%vph
        
#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)
        ! additional parameters for spatial profile
        ! e1, b1 from Lortentz transform
        current%e1 = current%e0 * (1-current%bf)
#ifdef APT_VACUUM_GAUSS
        current%b1 = current%b0 * (1-current%bf)        
#endif        
        current%rayl = 0.5*(1-current%bf)*current%w0**2*current%omega/c
        ! xf adjusted so focus will be at input xf at time t0
        current%xf = (1-current%bf)*current%xf - c*current%bf*current%t0 &
             + 0.5*c*current%bf*dt + current%bf*x_min
#endif
        
        ! t0 correction for evaluation time: t0 -> t0 - 0.5*dt
        ! additional correction so peak is at correct plane
        current%t0 = current%t0 - 0.5*dt - x_min/current%vg
        current%ph0 = current%ph0 + 0.5*current%omega*dt + current%omega*x_min/current%vph

#ifdef APT_PLASMA
        ! redefinition of wp2norm as amplitude of analytic current, for convenience
        current%wp2norm = -current%omega * (cvph2-1-(2*cvph2**2-(cvph2+1)*f2)*rfac/6.0) * current%e0 * epsilon0
#endif
        
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
    INTEGER :: n,i,err
    ! copy for convenience
    REAL(num) :: tenv_norm, itau_ivg, omega_ivph, pht
    REAL(num) :: e0_amp, b0_amp
    INTEGER :: sg
#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)    
    REAL(num) :: iw02, rayl, irayl, irayl2
    REAL(num) :: xf, yf, e1_amp
    REAL(num) :: amp, phase
    INTEGER :: j
#ifdef APT_VACUUM_GAUSS
    REAL(num) :: zf, b1_amp
    INTEGER :: k
#endif
#endif
#ifdef APT_PLASMA
    ! extra variables, used for finite duration correction
    REAL(num) :: env, phdt, adt1, adt2
#endif
    
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
          tenv_norm = (time-current%t0)/current%tau
          itau_ivg = 1.0_num/(current%tau*current%vg)
          sg = 2*current%sg
          omega_ivph = current%omega/current%vph
          pht = current%ph0 + current%omega*time
          e0_amp = current%e0
          b0_amp = current%b0
#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)    
          iw02 = 1.0_num/(current%w0*current%w0)
          rayl = current%rayl
          irayl = 1.0_num/current%rayl
          irayl2 = 1.0_num/(current%rayl*current%rayl)
          xf = current%xf + c*current%bf*time
          yf = current%yf
          e1_amp = current%e1
#ifdef APT_VACUUM_GAUSS
          zf = current%zf
          b1_amp = current%b1
#endif
#endif
          
#ifdef APT_PLASMA
          ! finite duration correction for analytic current case
          phdt = -((current%vg-current%vph)/current%vg)*sg/(current%omega*current%tau)
          adt1 = -0.5*((current%vg**2-current%vph**2)/current%vg**2)*sg**2/(current%omega*current%tau)**2
          adt2 = ((current%vg-current%vph)/current%vg)*sg*(sg-1)/(current%omega*current%tau)**2
#endif         
          
          ! ex (two parts) -- 3D or 2D gaussian only
#if defined(APT_VACUUM_GAUSS)
          DO k = 1-ng,nz+ng
            DO j = 1-ng,ny+ng
              DO i = 1-ng,nx+ng
                amp = e1_amp * irayl * (yb(j)-yf) &
                     * EXP( -iw02*((yb(j)-yf)*(yb(j)-yf)+(zb(k)-zf)*(zb(k)-zf))/(1+irayl2*(xf-x(i))**2) &
                     -(tenv_norm-itau_ivg*x(i))**sg ) / ((1+irayl2*(xf-x(i))**2)*SQRT(1+irayl2*(xf-x(i))**2))
                phase = pht-omega_ivph*x(i) - ATAN2(xf-x(i),rayl) &
                     + iw02*irayl*((yb(j)-yf)*(yb(j)-yf)+(zb(k)-zf)*(zb(k)-zf))*(xf-x(i))/(1+irayl2*(xf-x(i))**2)
                ex_total(i,j,k) = ex_total(i,j,k) + amp*( irayl*(xf-x(i))*SIN(phase) + COS(phase) )
              END DO
            END DO
          END DO
#elif defined(APT_VACUUM_GAUSS_2D)
          DO j = 1-ng,ny+ng
            DO i = 1-ng,nx+ng
              amp = e1_amp * irayl * (yb(j)-yf) &
                   * EXP( -iw02*(yb(j)-yf)*(yb(j)-yf)/(1+irayl2*(xf-x(i))**2) &
                   -(tenv_norm-itau_ivg*x(i))**sg ) / (1+irayl2*(xf-x(i))**2)**1.25
              phase = pht-omega_ivph*x(i) - 0.5*ATAN2(xf-x(i),rayl) &
                   + iw02*irayl*(yb(j)-yf)*(yb(j)-yf)*(xf-x(i))/(1+irayl2*(xf-x(i))**2)
              ex_total(i,j,:) = ex_total(i,j,:) + amp*( irayl*(xf-x(i))*SIN(phase) + COS(phase) )
            END DO
          END DO
#endif
          ! ey -- all cases
#if defined(APT_VACUUM_GAUSS)
          DO k = 1-ng,nz+ng
            DO j = 1-ng,ny+ng
              DO i = 1-ng,nx+ng
                ey_total(i,j,k) = ey_total(i,j,k) + e0_amp * &
                     EXP( -iw02*((y(j)-yf)*(y(j)-yf)+(zb(k)-zf)*(zb(k)-zf))/(1+irayl2*(xf-xb(i))**2) &
                     -(tenv_norm-itau_ivg*xb(i))**sg ) / SQRT(1+irayl2*(xf-xb(i))**2) &
                     *SIN( pht-omega_ivph*xb(i) - ATAN2(xf-xb(i),rayl) &
                     + iw02*irayl*((y(j)-yf)*(y(j)-yf)+(zb(k)-zf)*(zb(k)-zf))*(xf-xb(i))/(1+irayl2*(xf-xb(i))**2) )
              END DO
            END DO
          END DO
#elif defined(APT_VACUUM_GAUSS_2D)
          DO j = 1-ng,ny+ng
            DO i = 1-ng,nx+ng
              ey_total(i,j,:) = ey_total(i,j,:) + e0_amp * &
                   EXP( -iw02*(y(j)-yf)*(y(j)-yf)/(1+irayl2*(xf-xb(i))**2) &
                   -(tenv_norm-itau_ivg*xb(i))**sg ) / SQRT(SQRT(1+irayl2*(xf-xb(i))**2)) &
                   *SIN( pht-omega_ivph*xb(i) - 0.5*ATAN2(xf-xb(i),rayl) &
                   + iw02*irayl*(y(j)-yf)*(y(j)-yf)*(xf-xb(i))/(1+irayl2*(xf-xb(i))**2) )
            END DO
          END DO
#else
          DO i = 1-ng,nx+ng
            ey_total(i,:,:) = ey_total(i,:,:) + e0_amp * &
                 EXP( -(tenv_norm-itau_ivg*xb(i))**sg )*SIN( pht-omega_ivph*xb(i) )
          END DO          
#endif
          ! bx (two parts) -- 3D gaussian only
#if defined(APT_VACUUM_GAUSS)
          DO k = 1-ng,nz+ng
            DO j = 1-ng,ny+ng
              DO i = 1-ng,nx+ng
                amp = b1_amp * irayl * (z(k)-zf) &
                     * EXP( -iw02*((y(j)-yf)*(y(j)-yf)+(z(k)-zf)*(z(k)-zf))/(1+irayl2*(xf-xb(i))**2) &
                     -(tenv_norm-itau_ivg*xb(i))**sg ) / ((1+irayl2*(xf-xb(i))**2)*SQRT(1+irayl2*(xf-xb(i))**2))
                phase = pht-omega_ivph*xb(i) - ATAN2(xf-xb(i),rayl) &
                     + iw02*irayl*((y(j)-yf)*(y(j)-yf)+(z(k)-zf)*(z(k)-zf))*(xf-xb(i))/(1+irayl2*(xf-xb(i))**2)
                bx_total(i,j,k) = bx_total(i,j,k) + amp*( irayl*(xf-xb(i))*SIN(phase) + COS(phase) )
              END DO
            END DO
          END DO
#endif
          ! bz -- all cases
#if defined(APT_VACUUM_GAUSS)
          DO k = 1-ng,nz+ng
            DO j = 1-ng,ny+ng
              DO i = 1-ng,nx+ng
                bz_total(i,j,k) = bz_total(i,j,k) + b0_amp * &
                     EXP( -iw02*((y(j)-yf)*(y(j)-yf)+(zb(k)-zf)*(zb(k)-zf))/(1+irayl2*(xf-x(i))**2) &
                     -(tenv_norm-itau_ivg*x(i))**sg ) / SQRT(1+irayl2*(xf-x(i))**2) &
                     *SIN( pht-omega_ivph*x(i) - ATAN2(xf-x(i),rayl) &
                     + iw02*irayl*((y(j)-yf)*(y(j)-yf)+(zb(k)-zf)*(zb(k)-zf))*(xf-x(i))/(1+irayl2*(xf-x(i))**2) )
              END DO
            END DO
          END DO
#elif defined(APT_VACUUM_GAUSS_2D)
          DO j = 1-ng,ny+ng
            DO i = 1-ng,nx+ng
              bz_total(i,j,:) = bz_total(i,j,:) + b0_amp * &
                   EXP( -iw02*(y(j)-yf)*(y(j)-yf)/(1+irayl2*(xf-x(i))**2) &
                   -(tenv_norm-itau_ivg*x(i))**sg ) / SQRT(SQRT(1+irayl2*(xf-x(i))**2)) &
                   *SIN( pht-omega_ivph*x(i) - 0.5*ATAN2(xf-x(i),rayl) &
                   + iw02*irayl*(y(j)-yf)*(y(j)-yf)*(xf-x(i))/(1+irayl2*(xf-x(i))**2) )
            END DO
          END DO
#else
          DO i = 1-ng,nx+ng
#ifdef APT_VACUUM             
            bz_total(i,:,:) = bz_total(i,:,:) + b0_amp * &
                 EXP( -(tenv_norm-itau_ivg*x(i))**sg )*SIN( pht-omega_ivph*x(i) )
#else
            ! plasma case, needs finite duration correction
            env = tenv_norm-itau_ivg*x(i)
            bz_total(i,:,:) = bz_total(i,:,:) + b0_amp * &
                 EXP( -env**sg )*SIN( pht-omega_ivph*x(i) + phdt*env**(sg-1) ) &
                 *( 1 + adt1*env**(2*sg-2) + adt2*env**(sg-2) )
#endif            
          END DO
#endif
        END SELECT
        current => current%next
      END DO

  END SUBROUTINE analytic_pulse_update_arrays

  ! subroutines for use with analytic current
#ifdef APT_PLASMA
  SUBROUTINE analytic_pulse_update_j
    
    CALL analytic_pulse_update_j_analytic

    jx_diff = jx_diff + jx
    jy_diff = jy_diff + jy
    jz_diff = jz_diff + jz
    
  END SUBROUTINE analytic_pulse_update_j
  
  SUBROUTINE analytic_pulse_update_j_analytic

    TYPE(analytic_pulse_block), POINTER :: current
    TYPE(parameter_pack) :: parameters
    INTEGER :: n,i,err
    ! copy for convenience
    REAL(num) :: tenv_norm, itau_ivg, omega_ivph, pht
    REAL(num) :: j0_amp
    REAL(num) :: e0_amp, b0_amp, env
    INTEGER :: sg
    ! finite duration correction
    REAL(num) :: phdt, adt1, adt2

    ! important: j_diff = j_plasma - j_analytic
    ! this subroutine gives -j_analytic
    jx_diff = 0
    jy_diff = 0
    jz_diff = 0

    current => analytic_pulses
    
    DO WHILE(ASSOCIATED(current))
      ! evaluate the temporal evolution of the analytic_current
      ! loop over position
      ! note: time must be at n+1/2
      SELECT CASE(current%boundary)
          
      CASE(c_bd_x_min)         
        tenv_norm = (time-dt/2.0_num-current%t0)/current%tau
        itau_ivg = 1.0_num/(current%tau*current%vg)
        sg = 2*current%sg
        omega_ivph = current%omega/current%vph
        pht = current%ph0 + current%omega*(time-dt/2.0_num)
        j0_amp = current%wp2norm

        ! finite duration correction for current envelope
        phdt = -sg/(current%omega*current%tau)
        adt1 = -0.5*sg**2/(current%omega*current%tau)**2
        adt2 = sg*(sg-1)/(current%omega*current%tau)**2
         
        ! jy
        DO i = 1-jng,nx+jng
          env = tenv_norm-itau_ivg*xb(i)
          jy_diff(i,:,:) = jy_diff(i,:,:) + j0_amp * &
               EXP( -env**sg )*COS( pht-omega_ivph*xb(i) + phdt*env**(sg-1) ) &
               *( 1 + adt1*env**(2*sg-2) + adt2*env**(sg-2) )
        END DO
      END SELECT
      current => current%next
    END DO
    
  END SUBROUTINE analytic_pulse_update_j_analytic

  FUNCTION analytic_pulse_drift(part_pos,part_q,part_m)
    ! CAUTION: this feature is experimental and not extensively tested
    ! deterministic particle momentum for moving window or initialization
    ! note: requires a0 < 1 at particle location
    ! probably only good for a0 < 0.25 or so, since it ignores nonlinear response
    REAL(num), DIMENSION(c_ndirs) :: analytic_pulse_drift
    REAL(num), DIMENSION(c_ndims), INTENT(IN) :: part_pos
    REAL(num), INTENT(IN) :: part_q, part_m
    TYPE(analytic_pulse_block), POINTER :: current
    REAL(num), DIMENSION(c_ndirs) :: drift
    INTEGER :: n
    ! copy for convenience
    REAL(num) :: tenv_norm, itau_ivg, omega_ivph, pht, p_amp
    REAL(num) :: e0_amp, b0_amp, env, dmag
    INTEGER :: sg
    ! finite duration correction
    REAL(num) :: phdt, adt1, adt2, phdtc1, phdtc2
    
    drift = 0

    current => analytic_pulses
    
    DO WHILE(ASSOCIATED(current))
      SELECT CASE(current%boundary)          
      CASE(c_bd_x_min)         
        tenv_norm = (time-dt/2.0_num-current%t0)/current%tau
        itau_ivg = 1.0_num/(current%tau*current%vg)
        sg = 2*current%sg
        omega_ivph = current%omega/current%vph
        pht = current%ph0 + current%omega*(time-dt/2.0_num)
        ! - sign because part_q is -q0 for electron
        p_amp = -current%a0 * part_q/q0 * m0/part_m

        ! finite duration correction for current envelope
        phdt = -sg/(current%omega*current%tau)
        adt1 = -0.5*sg**2/(current%omega*current%tau)**2
        adt2 = sg*(sg-1)/(current%omega*current%tau)**2

        ! drift py
        env = tenv_norm-itau_ivg*part_pos(1)
        drift(2) = drift(2) + p_amp * &
             EXP( -env**sg )*COS( pht-omega_ivph*part_pos(1) + phdt*env**(sg-1) ) &
             *( 1 + adt1*env**(2*sg-2) + adt2*env**(sg-2) )

        current => current%next
      END SELECT
    END DO
    
    ! warning: will break horribly if drift > 1
    ! set max value just in case, gives p/mc ~ 2
    dmag = SQRT(SUM(drift*drift))
    drift = (drift/dmag)*MIN(dmag, 0.9_num)
    analytic_pulse_drift = part_m * c * drift/SQRT(1-SUM(drift*drift))
    
  END FUNCTION analytic_pulse_drift
#endif  
  
#endif
  
END MODULE analytic_pulse
