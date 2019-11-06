SUBROUTINE Subtract_CVs(state, sub)

! Subtract a control variable

USE DefConsTypes, ONLY :     &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT)   :: state
TYPE(CV_type),  INTENT(IN)      :: sub

state % v1(0:nlongs+1,0:nlevs+1) = state % v1(0:nlongs+1,0:nlevs+1) - &
                                   sub % v1(0:nlongs+1,0:nlevs+1)
state % v2(0:nlongs+1,0:nlevs+1) = state % v2(0:nlongs+1,0:nlevs+1) - &
                                   sub % v2(0:nlongs+1,0:nlevs+1)
state % v3(0:nlongs+1,0:nlevs+1) = state % v3(0:nlongs+1,0:nlevs+1) - &
                                   sub % v3(0:nlongs+1,0:nlevs+1)
state % v4(0:nlongs+1,0:nlevs+1) = state % v4(0:nlongs+1,0:nlevs+1) - &
                                   sub % v4(0:nlongs+1,0:nlevs+1)
state % v5(0:nlongs+1,0:nlevs+1) = state % v5(0:nlongs+1,0:nlevs+1) - &
                                   sub % v5(0:nlongs+1,0:nlevs+1)
state % v6(0:nlongs+1,0:nlevs+1) = state % v6(0:nlongs+1,0:nlevs+1) - &
                                   sub % v6(0:nlongs+1,0:nlevs+1)

END SUBROUTINE Subtract_CVs
