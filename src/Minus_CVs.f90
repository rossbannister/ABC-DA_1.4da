SUBROUTINE Minus_CVs(state)

! Replace the argument with minus its argument

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT)   :: state

REAL(ZREAL8)                    :: m1 = -1.0

state % v1(0:nlongs+1,0:nlevs+1) = m1 * state % v1(0:nlongs+1,0:nlevs+1)
state % v2(0:nlongs+1,0:nlevs+1) = m1 * state % v2(0:nlongs+1,0:nlevs+1)
state % v3(0:nlongs+1,0:nlevs+1) = m1 * state % v3(0:nlongs+1,0:nlevs+1)
state % v4(0:nlongs+1,0:nlevs+1) = m1 * state % v4(0:nlongs+1,0:nlevs+1)
state % v5(0:nlongs+1,0:nlevs+1) = m1 * state % v5(0:nlongs+1,0:nlevs+1)
state % v6(0:nlongs+1,0:nlevs+1) = m1 * state % v6(0:nlongs+1,0:nlevs+1)

END SUBROUTINE Minus_CVs
