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

MODULE deck_analytic_pulse_block

  USE strings_advanced
  USE analytic_pulse
  USE utilities

  IMPLICIT NONE
  SAVE

  PRIVATE
  PUBLIC :: analytic_pulse_deck_initialise, analytic_pulse_deck_finalise
  PUBLIC :: analytic_pulse_block_start, analytic_pulse_block_end
  PUBLIC :: analytic_pulse_block_handle_element, analytic_pulse_block_check

#ifdef APT_VACUUM
  TYPE(analytic_pulse_block), POINTER :: working_analytic_pulse
  LOGICAL :: boundary_set = .FALSE.
  INTEGER :: boundary
#endif

CONTAINS

  SUBROUTINE analytic_pulse_deck_initialise

#ifdef APT_VACUUM
    n_analytic_pulses(:) = 0
#endif

  END SUBROUTINE analytic_pulse_deck_initialise



  SUBROUTINE analytic_pulse_deck_finalise

#ifdef APT_VACUUM
    INTEGER :: io, iu
    ! check that user has only requested x_min boundary
    IF (SUM(n_analytic_pulses)-n_analytic_pulses(1) > 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Only analytic pulses launched from the x_min boundary are currently implemented'
          WRITE(io,*) 'Initialization from other boundaries will be ignored'
        END DO
    END IF
#endif

  END SUBROUTINE analytic_pulse_deck_finalise



  SUBROUTINE analytic_pulse_block_start

#ifdef APT_VACUUM
    IF (deck_state == c_ds_first) RETURN
    ! if any analytic pulse is defined
    use_analytic_pulses = .TRUE.
    ALLOCATE(working_analytic_pulse)
#endif

  END SUBROUTINE analytic_pulse_block_start



  SUBROUTINE analytic_pulse_block_end

#ifdef APT_VACUUM
    IF (deck_state == c_ds_first) RETURN

    CALL attach_analytic_pulse(working_analytic_pulse)
    boundary_set = .FALSE.
#endif

  END SUBROUTINE analytic_pulse_block_end



  FUNCTION analytic_pulse_block_handle_element(element, value) RESULT(errcode)

    CHARACTER(*), INTENT(IN) :: element, value
    INTEGER :: errcode
#ifdef APT_VACUUM
    REAL(num) :: dummy
    INTEGER :: io, iu
#endif

    errcode = c_err_none
#ifdef APT_VACUUM
    IF (deck_state == c_ds_first) RETURN
    IF (element == blank .OR. value == blank) RETURN

    IF (str_cmp(element, 'boundary')) THEN
      ! If the bounday has already been set, simply ignore further calls to it
      IF (boundary_set) RETURN
      boundary = as_boundary_print(value, element, errcode)
      boundary_set = .TRUE.
      CALL init_analytic_pulse(boundary, working_analytic_pulse)
      RETURN
    END IF

    IF (.NOT. boundary_set) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Cannot set analytic_pulse properties before boundary is set'
        END DO
        CALL abort_code(c_err_required_element_not_set)
      END IF
      extended_error_string = 'boundary'
      errcode = c_err_required_element_not_set
      RETURN
    END IF

    IF (str_cmp(element, 'omega') .OR. str_cmp(element, 'freq') &
         .OR. (str_cmp(element, 'frequency'))) THEN
      IF (rank == 0 .AND. str_cmp(element, 'freq')) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** WARNING ***'
          WRITE(io,*) 'Input deck line number ', TRIM(deck_line_number)
          WRITE(io,*) 'Element "freq" in the block "laser" is deprecated.'
          WRITE(io,*) 'Please use the element name "omega" instead.'
        END DO
      END IF
      working_analytic_pulse%omega = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'lambda')) THEN
      working_analytic_pulse%omega = 2*pi*c/as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'a0')) THEN
      working_analytic_pulse%a0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'tau')) THEN
      working_analytic_pulse%tau = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 't0')) THEN
      working_analytic_pulse%t0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'phase') .OR. str_cmp(element,'ph0')) THEN
      working_analytic_pulse%ph0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'sg')) THEN
      working_analytic_pulse%sg = as_integer_print(value, element, errcode)
      RETURN
    END IF

#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)
    ! parameters shared by 2D and 3D gaussians
    IF (str_cmp(element, 'w0')) THEN
      working_analytic_pulse%w0 = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'bf')) THEN
      working_analytic_pulse%bf = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'xf')) THEN
      working_analytic_pulse%xf = as_real_print(value, element, errcode)
      RETURN
    END IF

    IF (str_cmp(element, 'yf')) THEN
      working_analytic_pulse%yf = as_real_print(value, element, errcode)
      RETURN
    END IF

#ifdef APT_VACUUM_GAUSS
    ! only applicable to 3D gaussian
    IF (str_cmp(element, 'zf')) THEN
      working_analytic_pulse%zf = as_real_print(value, element, errcode)
      RETURN
    END IF
#endif
#endif

    errcode = c_err_unknown_element
#endif

  END FUNCTION analytic_pulse_block_handle_element



  FUNCTION analytic_pulse_block_check() RESULT(errcode)

    INTEGER :: errcode
#ifdef APT_VACUUM
    TYPE(analytic_pulse_block), POINTER :: current
    INTEGER :: error, io, iu
#endif

    errcode = c_err_none

#ifdef APT_VACUUM
    error = 0
    current => analytic_pulses
    DO WHILE(ASSOCIATED(current))
      IF (current%omega < 0.0_num) error = IOR(error, 1)
      IF (current%a0 < 0.0_num) error = IOR(error, 2)
      IF (current%tau < 0.0_num) error = IOR(error, 4)
#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)
      IF (current%w0 < 0.0_num) error = IOR(error, 8)
#endif
      current => current%next
    END DO

    IF (IAND(error, 1) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "lambda" or "omega" for every analytic pulse.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 2) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define an "a0" for every analytic pulse.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

    IF (IAND(error, 4) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "tau" for every analytic pulse.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF

#if defined(APT_VACUUM_GAUSS) || defined(APT_VACUUM_GAUSS_2D)
    IF (IAND(error, 8) /= 0) THEN
      IF (rank == 0) THEN
        DO iu = 1, nio_units ! Print to stdout and to file
          io = io_units(iu)
          WRITE(io,*) '*** ERROR ***'
          WRITE(io,*) 'Must define a "w0" for every analytic pulse.'
        END DO
      END IF
      errcode = c_err_missing_elements
    END IF
#endif

#endif

  END FUNCTION analytic_pulse_block_check

END MODULE deck_analytic_pulse_block
