SUBROUTINE Magnitude_rms (dataset, nxpoints, nypoints, rms)

! Code to calculate the RMS of a field

USE DefConsTypes, ONLY :  &
  ZREAL8

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: dataset(1:nxpoints,1:nypoints)
INTEGER,      INTENT(IN)    :: nxpoints, nypoints
REAL(ZREAL8), INTENT(INOUT) :: rms

INTEGER                     :: xx, yy
REAL(ZREAL8)                :: nxpoints_r, nypoints_r, ms, field


nxpoints_r = REAL(nxpoints)
nypoints_r = REAL(nypoints)

ms = 0.0
DO yy = 1, nypoints
  DO xx = 1, nxpoints
    field = dataset(xx,yy)
    ms    = ms + field * field
  END DO
END DO

ms  = ms / (nxpoints_r * nypoints_r)

rms = SQRT(ms)


END SUBROUTINE Magnitude_rms
