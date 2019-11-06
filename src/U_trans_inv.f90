SUBROUTINE U_trans_inv (LS, ControlVars, ModelVars, CVTdata, dims)

! Code to perform the inverse CVT: ControlVars = U_trans^-1 ModelVars

USE DefConsTypes, ONLY :  &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(ABC_type), INTENT(IN)    :: LS
TYPE(CV_type),  INTENT(INOUT) :: ControlVars
TYPE(ABC_type), INTENT(IN)    :: ModelVars
TYPE(CVT_type), INTENT(IN)    :: CVTdata
TYPE(dims_type),INTENT(IN)    :: dims

TYPE(CV_type)                 :: Interim1, Interim2


CALL U_p_inv (LS,                                                    &
              Interim2,                                              &
              ModelVars,                                             &
              CVTdata % CVT_order,                                   &
              CVTdata % CVT_param_opt_gb,                            &
              CVTdata % CVT_param_opt_hb,                            &
              CVTdata % CVT_param_opt_ab,                            &
              CVTdata % CVT_param_opt_reg,                           &
              CVTdata % Regression(1:nlevs,1:nlevs),                 &
              dims, .FALSE., '')

CALL U_stddev_inv (Interim2,                                         &
                   CVTdata)


! Determine the transform order
SELECT CASE (CVTdata % CVT_order)

CASE(1) ! The original Met Office transform order
        ! ---------------------------------------
  CALL U_v_inv (Interim1,                                            &
                Interim2,                                            &
                CVTdata)

  CALL U_h_inv (ControlVars,                                         &
                Interim1,                                            &
                CVTdata,                                             &
                .FALSE.)   ! (not just FFT)


CASE(2) ! The reversed horiz/vert transform order
        ! ---------------------------------------
  CALL U_h_inv (Interim1,                                            &
                Interim2,                                            &
                CVTdata,                                             &
                .FALSE.)   ! (not just FFT)

  CALL U_v_inv (ControlVars,                                         &
                Interim1,                                            &
                CVTdata)

END SELECT


END SUBROUTINE U_trans_inv
