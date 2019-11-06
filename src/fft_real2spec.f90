SUBROUTINE fft_real2spec (fieldr, fields)

! Code to do a balanced fft from real to spectral space

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  fft_worklen_x,          &
  fft_init_x,             &
  fft_wsave_x,            &
  fft_work_x,             &
  nlongs,                 &
  nlevs,                  &
  sr_nlongs,              &
  half_sr_nlongs

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: fieldr(1:nlongs,1:nlevs)
REAL(ZREAL8), INTENT(INOUT) :: fields(1:nlongs,1:nlevs)

INTEGER                     :: z, ierr

! Initialize the FFT in the x-direction if not already done
IF (.NOT.fft_init_x) THEN
  CALL rfft1i (nlongs, fft_wsave_x, fft_worklen_x, ierr)
  fft_init_x = .TRUE.
END IF


fields(1:nlongs,1:nlevs) = fieldr(1:nlongs,1:nlevs)
DO z = 1, nlevs
  ! Do fft in x for this z value
  CALL rfft1f (nlongs, 1,                                           &
               fields(1:nlongs,z), nlongs,                          &
               fft_wsave_x, fft_worklen_x, fft_work_x, nlongs, ierr)
  IF (ierr /= 0) THEN
    PRINT *, 'Error with rfft1f inside fft_real2spec'
    STOP
  END IF
END DO

!! Balance the FFT, and remove the 2 that appears in the library routine that we do not want
!fields(1,1:nlevs)          = fields(1,1:nlevs) * sr_nlongs
!fields(2:nlongs-1,1:nlevs) = fields(2:nlongs-1,1:nlevs) * half_sr_nlongs
!fields(nlongs,1:nlevs)     = fields(nlongs,1:nlevs) * sr_nlongs

! Balance the FFT
fields(1:nlongs,1:nlevs) = fields(1:nlongs,1:nlevs) * sr_nlongs

END SUBROUTINE fft_real2spec
