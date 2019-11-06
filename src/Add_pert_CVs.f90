SUBROUTINE Add_pert_CVs(state_out, state1_in, state2_in, weight)

! Computes the following
!   state_out = state1_in + weight * state2_in

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT)   :: state_out
TYPE(CV_type),  INTENT(IN)      :: state1_in
TYPE(CV_type),  INTENT(IN)      :: state2_in
REAL(ZREAL8),   INTENT(IN)      :: weight

state_out % v1(0:nlongs+1,0:nlevs+1) = state1_in % v1(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v1(0:nlongs+1,0:nlevs+1)

state_out % v2(0:nlongs+1,0:nlevs+1) = state1_in % v2(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v2(0:nlongs+1,0:nlevs+1)

state_out % v3(0:nlongs+1,0:nlevs+1) = state1_in % v3(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v3(0:nlongs+1,0:nlevs+1)

state_out % v4(0:nlongs+1,0:nlevs+1) = state1_in % v4(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v4(0:nlongs+1,0:nlevs+1)

state_out % v5(0:nlongs+1,0:nlevs+1) = state1_in % v5(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v5(0:nlongs+1,0:nlevs+1)

state_out % v6(0:nlongs+1,0:nlevs+1) = state1_in % v6(0:nlongs+1,0:nlevs+1) + &
                                       weight * state2_in % v6(0:nlongs+1,0:nlevs+1)

END SUBROUTINE Add_pert_CVs
