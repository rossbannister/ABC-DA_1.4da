SUBROUTINE U_trans (LS, ControlVars, ModelVars, CVTdata, dims)

! Code to perform the CVT: ModelVars = U ControlVars

USE DefConsTypes, ONLY :  &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs


IMPLICIT NONE

INCLUDE "Boundaries_CV.interface"


TYPE(ABC_type), INTENT(IN)    :: LS
TYPE(CV_type),  INTENT(IN)    :: ControlVars
TYPE(ABC_type), INTENT(INOUT) :: ModelVars
TYPE(CVT_type), INTENT(IN)    :: CVTdata
TYPE(dims_type),INTENT(IN)    :: dims

TYPE(CV_type)                 :: Interim1, Interim2


! Determine the transform order
SELECT CASE (CVTdata % CVT_order)

CASE(1) ! The original Met Office transform order
        ! ---------------------------------------
  CALL U_h (ControlVars,                                         &
            Interim1,                                            &
            CVTdata)

  CALL U_v (Interim1,                                            &
            Interim2,                                            &
            CVTdata)


CASE(2) ! The reversed horiz/vert transform order
        ! ---------------------------------------
  CALL U_v (ControlVars,                                         &
            Interim1,                                            &
            CVTdata)

  CALL U_h (Interim1,                                            &
            Interim2,                                            &
            CVTdata)

END SELECT

CALL U_stddev (Interim2,                                         &
               CVTdata)

CALL Boundaries_CV (Interim2)


CALL U_p (LS,                                                    &
          Interim2,                                              &
          ModelVars,                                             &
          CVTdata % CVT_order,                                   &
          CVTdata % CVT_param_opt_gb,                            &
          CVTdata % CVT_param_opt_hb,                            &
          CVTdata % CVT_param_opt_ab,                            &
          CVTdata % CVT_param_opt_reg,                           &
          CVTdata % Regression(1:nlevs,1:nlevs),                 &
          dims)

END SUBROUTINE U_trans
