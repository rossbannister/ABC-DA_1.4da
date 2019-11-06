SUBROUTINE Add_model_vars(state, increment, core_only)

! Increment a state

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(ABC_type), INTENT(INOUT)   :: state    ! Running total
TYPE(ABC_type), INTENT(IN)      :: increment
LOGICAL,        INTENT(IN)      :: core_only


  state % u(0:nlongs+1,0:nlevs+1)           = state % u(0:nlongs+1,0:nlevs+1) + &
                                              increment % u(0:nlongs+1,0:nlevs+1)
  state % v(0:nlongs+1,0:nlevs+1)           = state % v(0:nlongs+1,0:nlevs+1) + &
                                              increment % v(0:nlongs+1,0:nlevs+1)
  state % w(0:nlongs+1,0:nlevs+1)           = state % w(0:nlongs+1,0:nlevs+1) + &
                                              increment % w(0:nlongs+1,0:nlevs+1)
  state % r(0:nlongs+1,0:nlevs+1)           = state % r(0:nlongs+1,0:nlevs+1) + &
                                              increment % r(0:nlongs+1,0:nlevs+1)
  state % b(0:nlongs+1,0:nlevs+1)           = state % b(0:nlongs+1,0:nlevs+1) + &
                                              increment % b(0:nlongs+1,0:nlevs+1)

  IF (.NOT.core_only) THEN
    state % rho(0:nlongs+1,0:nlevs+1)         = state % rho(0:nlongs+1,0:nlevs+1) + &
                                                increment % rho(0:nlongs+1,0:nlevs+1)
    state % b_ef(0:nlongs+1,0:nlevs+1)        = state % b_ef(0:nlongs+1,0:nlevs+1) + &
                                                increment % b_ef(0:nlongs+1,0:nlevs+1)
    state % tracer(0:nlongs+1,0:nlevs+1)      = state % tracer(0:nlongs+1,0:nlevs+1) + &
                                                increment % tracer(0:nlongs+1,0:nlevs+1)
    state % hydro_imbal(0:nlongs+1,0:nlevs+1) = state % hydro_imbal(0:nlongs+1,0:nlevs+1) + &
                                                increment % hydro_imbal(0:nlongs+1,0:nlevs+1)
    state % geost_imbal(0:nlongs+1,0:nlevs+1) = state % geost_imbal(0:nlongs+1,0:nlevs+1) + &
                                                increment % geost_imbal(0:nlongs+1,0:nlevs+1)
    state % Kinetic_Energy                    = state % Kinetic_Energy + &
                                                increment % Kinetic_Energy
    state % Buoyant_Energy                    = state % Buoyant_Energy + &
                                                increment % Buoyant_Energy
    state % Elastic_Energy                    = state % Elastic_Energy + &
                                                increment % Elastic_Energy
    state % Total_Energy                      = state % Total_Energy + &
                                                increment % Total_Energy
  END IF

END SUBROUTINE Add_model_vars
