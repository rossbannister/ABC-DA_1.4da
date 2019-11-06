SUBROUTINE Initialise_model_vars(state, random_init)

! Initialise model variables

USE DefConsTypes, ONLY :     &
    ABC_type,                &
    nlongs, nlevs

IMPLICIT NONE

INCLUDE "Boundaries.interface"

TYPE(ABC_type), INTENT(INOUT)   :: state
LOGICAL,        INTENT(IN)      :: random_init

IF (random_init) THEN

  CALL RANDOM_NUMBER (state % u(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % v(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % w(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % r(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % b(1:nlongs,1:nlevs))
  CALL RANDOM_NUMBER (state % tracer(1:nlongs,1:nlevs))

  CALL Boundaries (state)

ELSE

  state % u(0:nlongs+1,0:nlevs+1)           = 0.0
  state % v(0:nlongs+1,0:nlevs+1)           = 0.0
  state % w(0:nlongs+1,0:nlevs+1)           = 0.0
  state % r(0:nlongs+1,0:nlevs+1)           = 0.0     ! density perturbation
  state % b(0:nlongs+1,0:nlevs+1)           = 0.0     ! buoyancy perturbation
  state % rho(0:nlongs+1,0:nlevs+1)         = 0.0     ! density full field
  state % b_ef(0:nlongs+1,0:nlevs+1)        = 0.0     ! Effective buoyancy
  state % tracer(0:nlongs+1,0:nlevs+1)      = 0.0
  state % hydro_imbal(0:nlongs+1,0:nlevs+1) = 0.0
  state % geost_imbal(0:nlongs+1,0:nlevs+1) = 0.0
  state % Kinetic_Energy                    = 0.0
  state % Buoyant_Energy                    = 0.0
  state % Elastic_Energy                    = 0.0
  state % Total_Energy                      = 0.0

END IF

END SUBROUTINE Initialise_model_vars
