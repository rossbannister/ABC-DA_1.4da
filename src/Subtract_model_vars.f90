SUBROUTINE Subtract_model_vars(state, sub, core_only)

! Subtract a state

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state
TYPE(ABC_type), INTENT(IN)      :: sub
LOGICAL,        INTENT(IN)      :: core_only


  state % u(0:nlongs+1,0:nlevs+1)           = state % u(0:nlongs+1,0:nlevs+1) - &
                                              sub % u(0:nlongs+1,0:nlevs+1)
  state % v(0:nlongs+1,0:nlevs+1)           = state % v(0:nlongs+1,0:nlevs+1) - &
                                              sub % v(0:nlongs+1,0:nlevs+1)
  state % w(0:nlongs+1,0:nlevs+1)           = state % w(0:nlongs+1,0:nlevs+1) - &
                                              sub % w(0:nlongs+1,0:nlevs+1)
  state % r(0:nlongs+1,0:nlevs+1)           = state % r(0:nlongs+1,0:nlevs+1) - &
                                              sub % r(0:nlongs+1,0:nlevs+1)
  state % b(0:nlongs+1,0:nlevs+1)           = state % b(0:nlongs+1,0:nlevs+1) - &
                                              sub % b(0:nlongs+1,0:nlevs+1)
  state % tracer(0:nlongs+1,0:nlevs+1)      = state % tracer(0:nlongs+1,0:nlevs+1) - &
                                              sub % tracer(0:nlongs+1,0:nlevs+1)

  IF (.NOT.core_only) THEN
    state % rho(0:nlongs+1,0:nlevs+1)         = state % rho(0:nlongs+1,0:nlevs+1) - &
                                                sub % rho(0:nlongs+1,0:nlevs+1)
    state % b_ef(0:nlongs+1,0:nlevs+1)        = state % b_ef(0:nlongs+1,0:nlevs+1) - &
                                                sub % b_ef(0:nlongs+1,0:nlevs+1)
    state % hydro_imbal(0:nlongs+1,0:nlevs+1) = state % hydro_imbal(0:nlongs+1,0:nlevs+1) - &
                                                sub % hydro_imbal(0:nlongs+1,0:nlevs+1)
    state % geost_imbal(0:nlongs+1,0:nlevs+1) = state % geost_imbal(0:nlongs+1,0:nlevs+1) - &
                                                sub % geost_imbal(0:nlongs+1,0:nlevs+1)
    state % Kinetic_Energy                    = state % Kinetic_Energy - &
                                                sub % Kinetic_Energy
    state % Buoyant_Energy                    = state % Buoyant_Energy - &
                                                sub % Buoyant_Energy
    state % Elastic_Energy                    = state % Elastic_Energy - &
                                                sub % Elastic_Energy
    state % Total_Energy                      = state % Total_Energy - &
                                                sub % Total_Energy
  END IF

END SUBROUTINE Subtract_model_vars
