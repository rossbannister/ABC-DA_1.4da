SUBROUTINE U_trans_adj (LS, ControlVars, ModelVars, CVTdata, dims)

! Code to perform the adjoint of the CVT: ControlVars = U* ModelVars

USE DefConsTypes, ONLY :  &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs


IMPLICIT NONE

INCLUDE "Boundaries_CV_adj.interface"


TYPE(ABC_type), INTENT(IN)    :: LS
TYPE(CV_type),  INTENT(INOUT) :: ControlVars
TYPE(ABC_type), INTENT(IN)    :: ModelVars
TYPE(CVT_type), INTENT(IN)    :: CVTdata
TYPE(dims_type),INTENT(IN)    :: dims

TYPE(CV_type)                 :: Interim1, Interim2



CALL Initialise_CVs(Interim2, .FALSE.)

CALL U_p_adj (LS,                                                    &
              Interim2,                                              &
              ModelVars,                                             &
              CVTdata % CVT_order,                                   &
              CVTdata % CVT_param_opt_gb,                            &
              CVTdata % CVT_param_opt_hb,                            &
              CVTdata % CVT_param_opt_ab,                            &
              CVTdata % CVT_param_opt_reg,                           &
              CVTdata % Regression(1:nlevs,1:nlevs),                 &
              dims)

CALL Boundaries_CV_adj (Interim2)


CALL U_stddev (Interim2,                                             &
               CVTdata)

! Determine the transform order
SELECT CASE (CVTdata % CVT_order)

CASE(1) ! The original Met Office transform order
        ! ---------------------------------------
  CALL U_v_adj (Interim1,                                            &
                Interim2,                                            &
                CVTdata)

  CALL U_h_adj (ControlVars,                                         &
                Interim1,                                            &
                CVTdata)


CASE(2) ! The reversed horiz/vert transform order
        ! ---------------------------------------
  CALL U_h_adj (Interim1,                                            &
                Interim2,                                            &
                CVTdata)

  CALL U_v_adj (ControlVars,                                         &
                Interim1,                                            &
                CVTdata)

END SELECT


END SUBROUTINE U_trans_adj
