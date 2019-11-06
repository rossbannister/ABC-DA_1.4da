SUBROUTINE Div_model_cons(state, div_cons, core_only)

! Divide a state by a constant

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state
REAL(ZREAL8),   INTENT(IN)      :: div_cons
LOGICAL,        INTENT(IN)      :: core_only

REAL(ZREAL8)                    :: mul_cons


  mul_cons = 1.0 / div_cons

  state % u(0:nlongs+1,0:nlevs+1)           = state % u(0:nlongs+1,0:nlevs+1) * mul_cons
  state % v(0:nlongs+1,0:nlevs+1)           = state % v(0:nlongs+1,0:nlevs+1) * mul_cons
  state % w(0:nlongs+1,0:nlevs+1)           = state % w(0:nlongs+1,0:nlevs+1) * mul_cons
  state % r(0:nlongs+1,0:nlevs+1)           = state % r(0:nlongs+1,0:nlevs+1) * mul_cons
  state % b(0:nlongs+1,0:nlevs+1)           = state % b(0:nlongs+1,0:nlevs+1) * mul_cons

  IF (.NOT.core_only) THEN
    state % rho(0:nlongs+1,0:nlevs+1)         = state % rho(0:nlongs+1,0:nlevs+1) * mul_cons
    state % b_ef(0:nlongs+1,0:nlevs+1)        = state % b_ef(0:nlongs+1,0:nlevs+1) * mul_cons
    state % tracer(0:nlongs+1,0:nlevs+1)      = state % tracer(0:nlongs+1,0:nlevs+1) * mul_cons
    state % hydro_imbal(0:nlongs+1,0:nlevs+1) = state % hydro_imbal(0:nlongs+1,0:nlevs+1) * mul_cons
    state % geost_imbal(0:nlongs+1,0:nlevs+1) = state % geost_imbal(0:nlongs+1,0:nlevs+1) * mul_cons
    state % Kinetic_Energy                    = state % Kinetic_Energy * mul_cons
    state % Buoyant_Energy                    = state % Buoyant_Energy * mul_cons
    state % Elastic_Energy                    = state % Elastic_Energy * mul_cons
    state % Total_Energy                      = state % Total_Energy * mul_cons
  END IF

END SUBROUTINE Div_model_cons
