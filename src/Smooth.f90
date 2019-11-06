SUBROUTINE Smooth (size, field, points)

!********************************************************
!* Subroutine to smooth the input array                 *
!*                                                      *
!* R. Bannister, vn1.4da, 29-07-19                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8

IMPLICIT NONE
INTEGER,        INTENT(IN)    :: size
REAL(ZREAL8),   INTENT(INOUT) :: field(1:size)
INTEGER,        INTENT(IN)    :: points

INTEGER                       :: l1, l2
INTEGER                       :: start, finish
REAL(ZREAL8)                  :: npoints, total
REAL(ZREAL8)                  :: smoothed(1:size)

! Do the smoothing
DO l1 = 1, size
  start        = l1 - points
  IF (start < 1) start = 1
  finish       = l1 + points
  IF (finish > size) finish = size
  npoints      = REAL(finish - start + 1)
  total        = SUM(field(start:finish))
  smoothed(l1) = total / npoints
END DO

! Overwrite the input
field(1:size) = smoothed(1:size)

END SUBROUTINE Smooth
