SUBROUTINE LinearBal_r_adj (psi, r_b)

! Code to compute the adjoint version of the balanced component of r

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  f, C, half

IMPLICIT NONE

REAL(ZREAL8),   INTENT(INOUT) :: psi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: r_b(1:nlongs,1:nlevs)

INTEGER                       :: x
REAL(ZREAL8)                  :: fac, d(1:nlevs)

!See Eq (18a) of model paper
fac = f * half / C
DO x = 1, nlongs
  d(1:nlevs)       = fac * r_b(x,1:nlevs)
  psi(x,1:nlevs)   = psi(x,1:nlevs) + d
  psi(x-1,1:nlevs) = psi(x-1,1:nlevs) + d
END DO


END SUBROUTINE LinearBal_r_adj
