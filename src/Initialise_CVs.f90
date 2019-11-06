SUBROUTINE Initialise_CVs(state, random_init)

! Initialise control variables

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CV_type,                 &
    nlongs, nlevs,           &
    unity

IMPLICIT NONE

INCLUDE "Boundaries_CV.interface"

TYPE(CV_type),  INTENT(INOUT)   :: state
LOGICAL,        INTENT(IN)      :: random_init

INTEGER                         :: x, z
REAL(ZREAL8)                    :: GAUSS

IF (random_init) THEN

  DO x = 1, nlongs
    DO z = 1, nlevs
      state % v1(x,z) = GAUSS (unity)
      state % v2(x,z) = GAUSS (unity)
      state % v3(x,z) = GAUSS (unity)
      state % v4(x,z) = GAUSS (unity)
      state % v5(x,z) = GAUSS (unity)
      state % v6(x,z) = GAUSS (unity)
    END DO
  END DO

  CALL Boundaries_CV (state)

ELSE

  state % v1(0:nlongs+1,0:nlevs+1) = 0.0
  state % v2(0:nlongs+1,0:nlevs+1) = 0.0
  state % v3(0:nlongs+1,0:nlevs+1) = 0.0
  state % v4(0:nlongs+1,0:nlevs+1) = 0.0
  state % v5(0:nlongs+1,0:nlevs+1) = 0.0
  state % v6(0:nlongs+1,0:nlevs+1) = 0.0

END IF

END SUBROUTINE Initialise_CVs
