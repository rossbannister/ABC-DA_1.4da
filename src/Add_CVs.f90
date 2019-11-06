SUBROUTINE Add_CVs(state, increment)

! Add control variables

USE DefConsTypes, ONLY :     &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT)   :: state
TYPE(CV_type),  INTENT(IN)      :: increment

state % v1(0:nlongs+1,0:nlevs+1) = state % v1(0:nlongs+1,0:nlevs+1) + &
                                   increment % v1(0:nlongs+1,0:nlevs+1)
state % v2(0:nlongs+1,0:nlevs+1) = state % v2(0:nlongs+1,0:nlevs+1) + &
                                   increment % v2(0:nlongs+1,0:nlevs+1)
state % v3(0:nlongs+1,0:nlevs+1) = state % v3(0:nlongs+1,0:nlevs+1) + &
                                   increment % v3(0:nlongs+1,0:nlevs+1)
state % v4(0:nlongs+1,0:nlevs+1) = state % v4(0:nlongs+1,0:nlevs+1) + &
                                   increment % v4(0:nlongs+1,0:nlevs+1)
state % v5(0:nlongs+1,0:nlevs+1) = state % v5(0:nlongs+1,0:nlevs+1) + &
                                   increment % v5(0:nlongs+1,0:nlevs+1)
state % v6(0:nlongs+1,0:nlevs+1) = state % v6(0:nlongs+1,0:nlevs+1) + &
                                   increment % v6(0:nlongs+1,0:nlevs+1)

END SUBROUTINE Add_CVs
