SUBROUTINE HydroBal_b_adj (r, b_b, dims)

! Code to compute the ajoint version of the hydrostatically balanced component of b

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  C

IMPLICIT NONE

REAL(ZREAL8),    INTENT(INOUT) :: r(1:nlongs,0:nlevs+1)
REAL(ZREAL8),    INTENT(IN)    :: b_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

INTEGER                        :: z
REAL(ZREAL8)                   :: recipdz
REAL(ZREAL8)                   :: d(1:nlongs)

!See Eq (19) of model paper
DO z = 1, nlevs
  recipdz = 1.0 / (dims % half_levs(z+1) - dims % half_levs(z))
  d(1:nlongs)     = C * recipdz * b_b(1:nlongs,z)
  r(1:nlongs,z+1) = r(1:nlongs,z+1) + d(1:nlongs)
  r(1:nlongs,z)   = r(1:nlongs,z) - d(1:nlongs)
END DO


END SUBROUTINE HydroBal_b_adj
