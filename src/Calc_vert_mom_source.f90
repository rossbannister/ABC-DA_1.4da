SUBROUTINE Calc_vert_mom_source (state, dims)

! Calculate source of vertical momentum in real space

USE DefConsTypes, ONLY :   &
    ABC_type,              &
    dims_type,             &
    ZREAL8,                &
    nlongs, nlevs, C

IMPLICIT NONE

! Declare parameters
!---------------------
TYPE(ABC_type), INTENT(INOUT)  :: state
TYPE(dims_type), INTENT(IN)    :: dims

! Declare variables
!---------------------
INTEGER                               :: x, z

! This has a very similiar form to the hydrostatic imbalance

DO x = 1, nlongs
  DO z = 1, nlevs
    state % vert_mom_source(x,z) = ( state % b(x,z) -                                  &
                                     C * (state % r(x,z+1) - state % r(x,z)) /         &
                                     (dims % half_levs(z+1) - dims % half_levs(z)) ) * &
                                   state % rho(x,z)
  ENDDO
ENDDO


END SUBROUTINE Calc_vert_mom_source
