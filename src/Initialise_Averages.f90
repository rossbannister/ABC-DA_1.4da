SUBROUTINE Initialise_Averages(state)

! Initialise averages variables

USE DefConsTypes, ONLY :     &
    Averages_type,           &
    nlongs, nlevs

IMPLICIT NONE

TYPE(Averages_type), INTENT(INOUT)   :: state

state % u_1(0:nlongs+1, 0:nlevs+1) = 0.0
state % u_2(0:nlongs+1, 0:nlevs+1) = 0.0
state % u_m(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_1(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_2(0:nlongs+1, 0:nlevs+1) = 0.0
state % w_m(0:nlongs+1, 0:nlevs+1) = 0.0

END SUBROUTINE Initialise_Averages
