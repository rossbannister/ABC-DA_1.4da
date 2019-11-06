SUBROUTINE Initialise_um_data (state)

! Initialise um-data variables

USE DefConsTypes, ONLY :   &
    UM_type,               &
    nlongs, nlevs

IMPLICIT NONE

TYPE(UM_type), INTENT(INOUT)   :: state

state % longs_u(1:nlongs)                  = 0.0
state % longs_v(1:nlongs)                  = 0.0
state % half_levs(1:nlevs+1)               = 0.0
state % full_levs(0:nlevs)                 = 0.0
state % u(1:nlongs,1:nlevs)                = 0.0
state % v(1:nlongs,1:nlevs)                = 0.0
state % w(1:nlongs,0:nlevs)                = 0.0
state % density(1:nlongs,1:nlevs)          = 0.0
state % theta(1:nlongs,1:nlevs)            = 0.0
state % exner_pressure(1:nlongs,1:nlevs+1) = 0.0
state % orog_height(1:nlongs)              = 0.0

END SUBROUTINE Initialise_um_data
