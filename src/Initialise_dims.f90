SUBROUTINE Initialise_dims (state)

! Initialise dimension variables

USE DefConsTypes, ONLY :   &
    dims_type,             &
    nlongs, nlevs

IMPLICIT NONE

TYPE(dims_type), INTENT(INOUT)   :: state

state % longs_u(0:nlongs+1)  = 0.0
state % longs_v(0:nlongs+1)  = 0.0
state % half_levs(0:nlevs+1) = 0.0
state % full_levs(0:nlevs+1) = 0.0
state % a1(1:nlevs)          = 0.0
state % b1(1:nlevs)          = 0.0
state % a2(1:nlevs)          = 0.0
state % b2(1:nlevs)          = 0.0

END SUBROUTINE Initialise_dims
