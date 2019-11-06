SUBROUTINE fft_spec2real (fields, fieldr)

! Code to do a balanced fft from spectral to real space

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  fft_worklen_x,          &
  fft_init_x,             &
  fft_wsave_x,            &
  fft_work_x,             &
  nlongs,                 &
  nlevs,                  &
  sr_nlongs

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: fields(1:nlongs,1:nlevs)
REAL(ZREAL8), INTENT(INOUT) :: fieldr(1:nlongs,1:nlevs)

INTEGER                     :: z, ierr

! Initialize the FFT in the x-direction if not already done
IF (.NOT.fft_init_x) THEN
  CALL rfft1i (nlongs, fft_wsave_x, fft_worklen_x, ierr)
  fft_init_x = .TRUE.
END IF

!! Undo the half that was done at the end of fft_real2spec
!fieldr(1,1:nlevs)          = fields(1,1:nlevs)
!fieldr(2:nlongs-1,1:nlevs) = 2. * fields(2:nlongs-1,1:nlevs)
!fieldr(nlongs,1:nlevs)     = fields(nlongs,1:nlevs)


fieldr(1:nlongs,1:nlevs) = fields(1:nlongs,1:nlevs)

DO z = 1, nlevs
  ! Do inverse fft in x for this z value
  CALL rfft1b (nlongs, 1,                                           &
               fieldr(1:nlongs,z), nlongs,                          &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  IF (ierr /= 0) THEN
    PRINT *, 'Error with rfft1b inside fft_spec2real'
    STOP
  END IF
END DO

fieldr(1:nlongs,1:nlevs) = fieldr(1:nlongs,1:nlevs) / sr_nlongs


END SUBROUTINE fft_spec2real
