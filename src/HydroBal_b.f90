SUBROUTINE HydroBal_b (r, b_b, dims)

! Code to compute the hydrostatically balanced component of b

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  C

IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)    :: r(1:nlongs,0:nlevs+1)
REAL(ZREAL8),    INTENT(INOUT) :: b_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

REAL(ZREAL8)                   :: recipdz
INTEGER                        :: z

!See Eq (19) of model paper
DO z = 1, nlevs
  recipdz = 1.0 / (dims % half_levs(z+1) - dims % half_levs(z))
  b_b(1:nlongs,z) = C * recipdz * (r(1:nlongs,z+1) - r(1:nlongs,z))
END DO


END SUBROUTINE HydroBal_b
