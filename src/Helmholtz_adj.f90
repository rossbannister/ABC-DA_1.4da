SUBROUTINE Helmholtz_adj (psi, chi, u, v)

! Code to perform the Helmholtz (recover u and v from psi and chi)

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  recipdx

IMPLICIT NONE

REAL(ZREAL8),   INTENT(INOUT) :: psi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(INOUT) :: chi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: u(1:nlongs,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: v(1:nlongs,1:nlevs)

INTEGER                       :: x

! See Eq (39) of the model paper
DO x = 1, nlongs
  chi(x+1,1:nlevs) = chi(x+1,1:nlevs) + u(x,1:nlevs) * recipdx
  chi(x,1:nlevs)   = chi(x,1:nlevs) - u(x,1:nlevs) * recipdx
  psi(x,1:nlevs)   = psi(x,1:nlevs) + v(x,1:nlevs) * recipdx
  psi(x-1,1:nlevs) = psi(x-1,1:nlevs) - v(x,1:nlevs) * recipdx
END DO


END SUBROUTINE Helmholtz_adj
