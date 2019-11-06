SUBROUTINE LinearBal_r (psi, r_b)

! Code to compute the balanced component of r

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  f, C, half

IMPLICIT NONE

REAL(ZREAL8),   INTENT(IN)    :: psi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(INOUT) :: r_b(1:nlongs,1:nlevs)

INTEGER                       :: x
REAL(ZREAL8)                  :: fac

!See Eq (18a) of model paper
fac = f * half / C
DO x = 1, nlongs
  r_b(x,1:nlevs) = fac * (psi(x,1:nlevs) + psi(x-1,1:nlevs))
END DO


END SUBROUTINE LinearBal_r
