SUBROUTINE Calc_horiz_div (state, Dims)

! Calculate the horizontal divergence in real space

USE DefConsTypes, ONLY :   &
    ABC_type,              &
    dims_type,             &
    ZREAL8,                &
    nlongs, nlevs, recipdx

IMPLICIT NONE

! Declare parameters
!---------------------
TYPE(ABC_type), INTENT(INOUT)  :: state
TYPE(dims_type), INTENT(IN)    :: Dims

! Declare variables
!---------------------
INTEGER                        :: x, z

DO x = 1, nlongs
  DO z = 1, nlevs
    state % horiz_div(x,z) = ( state % u(x,z) - state % u(x-1,z) ) * recipdx
  ENDDO
ENDDO


END SUBROUTINE Calc_horiz_div
