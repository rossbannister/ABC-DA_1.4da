SUBROUTINE Helmholtz (psi, chi, u, v)

! Code to perform the Helmholtz (recover u and v from psi and chi)

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  recipdx

IMPLICIT NONE

REAL(ZREAL8),   INTENT(IN)    :: psi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(IN)    :: chi(0:nlongs+1,1:nlevs)
REAL(ZREAL8),   INTENT(INOUT) :: u(1:nlongs,1:nlevs)
REAL(ZREAL8),   INTENT(INOUT) :: v(1:nlongs,1:nlevs)

INTEGER                       :: x

! See Eq (39) of the model paper
DO x = 1, nlongs
  u(x,1:nlevs) = (chi(x+1,1:nlevs) - chi(x,1:nlevs)) * recipdx
  v(x,1:nlevs) = (psi(x,1:nlevs) - psi(x-1,1:nlevs)) * recipdx
END DO


END SUBROUTINE Helmholtz
