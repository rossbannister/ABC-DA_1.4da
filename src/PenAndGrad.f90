SUBROUTINE PenAndGrad ( VarType, dasteps, steps, dims, t0,                       &
                        Obs, dchi0, dchi0zero, dchib0, dchib0zero,               &
                        CVTdata, LSfc, Jb, Jo, J, compute_grad, grad0,           &
                        outputdiags_flag, forceNLmodel, outputdir, outerloop, innerloop )

!*********************************************************************************
!*                                                                               *
!*  Compute penalty and gradient of the cost function                            *
!*  The gradient is in control space                                             *
!*                                                                               *
!*  INPUTS                                                                       *
!*  VarType                - 3=3DVar, 35=3D-FGAT, 4=4DVar                        *
!*  dasteps                - number of time states for the da (excluding 0)      *
!*  steps                  - number of model steps                               *
!*  dims                   - dimension data                                      *
!*  t0                     - time corresponding to time step 0 (seconds)         *
!*  Obs                    - pointer to the start of the observation linked list *
!*  dchi0                  - control variable at start of window                 *
!*  dchi0zero              - set if dchi0=0                                      *
!*  dchib0                 - background increment in control space               *
!*  dchib0zero             - set if dchib0=0                                     *
!*  CVTdata                - data that defines the CVT                           *
!*  INPUT/OUTPUT                                                                 *
!*  LSfc                   - Linearization state trajectory                      *
!*  OUTPUTS                                                                      *
!*  Jb                     - Background part of cost function                    *
!*  Jo                     - Observation part of cost function                   *
!*  J                      - Total cost function                                 *
!*  compute_grad           - Set to compute the gradient (otherwise only cost    *
!*                           function is found, leaving grad0 unchanged)         *
!*  grad0                  - Gradient of cost function in control space          *
!*  INPUTS                                                                       *
!*  outputdiags_flag       - Set to output various optional diagnostics          *
!*  forceNLmodel           - Set to force run of NL model in ref state           *
!*  outputdir              - Output directory                                    *
!*  outerloop              - outer loop index (for diagnostic file naming)       *
!*  innerloop              - inner loop index (for diagnostic file naming)       *
!*                                                                               *
!*   R. Bannister, 1.4da 28-02-2018                                              *
!*                                                                               *
!*********************************************************************************


USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  Obs_type,               &
  CV_type,                &
  CVT_type,               &
  dt, dt_da

IMPLICIT NONE

INCLUDE "InnerProduct_obsself.interface"
INCLUDE "InnerProdControlSpace.interface"
INCLUDE "ModelObservations.interface"
INCLUDE "ModelObservations_ZeroPert.interface"
INCLUDE "ModelObservations_linear.interface"
INCLUDE "ModelObservations_adj.interface"
INCLUDE "Write_Obs.interface"

! Parameters
!-----------
INTEGER,                 INTENT(IN)    :: Vartype
INTEGER,                 INTENT(IN)    :: dasteps
INTEGER,                 INTENT(IN)    :: steps
TYPE(dims_type),         INTENT(IN)    :: dims
INTEGER,                 INTENT(IN)    :: t0
TYPE(Obs_type), POINTER, INTENT(IN)    :: Obs
TYPE(CV_type),           INTENT(IN)    :: dchi0
LOGICAL,                 INTENT(IN)    :: dchi0zero
TYPE(CV_type),           INTENT(IN)    :: dchib0
LOGICAL,                 INTENT(IN)    :: dchib0zero
TYPE(CVT_type),          INTENT(IN)    :: CVTdata
TYPE(ABC_type),          INTENT(INOUT) :: LSfc(0:steps)
REAL(ZREAL8),            INTENT(OUT)   :: Jb
REAL(ZREAL8),            INTENT(OUT)   :: Jo
REAL(ZREAL8),            INTENT(OUT)   :: J
LOGICAL,                 INTENT(IN)    :: compute_grad
TYPE(CV_type),           INTENT(INOUT) :: grad0
LOGICAL,                 INTENT(IN)    :: outputdiags_flag
LOGICAL,                 INTENT(IN)    :: forceNLmodel
CHARACTER(LEN=100),      INTENT(IN)    :: outputdir
INTEGER,                 INTENT(IN)    :: innerloop
INTEGER,                 INTENT(IN)    :: outerloop

! Local variables
! ---------------
TYPE(ABC_type)                         :: deltax(0:dasteps)
TYPE(ABC_type)                         :: deltax_hat(0:dasteps)
INTEGER                                :: t, ModelStepsPerDAStep
COMPLEX(ZREAL8)                        :: Jb_complex
TYPE(CV_type)                          :: diffcv
CHARACTER(LEN=320)                     :: Diag_filename



ModelStepsPerDAStep = INT(dt_da/dt)

! ===================================
! --- FORWARD PART OF CALCULATION ---
! ===================================
IF (dchi0zero) THEN
  ! This means that we are at the start of a new outer loop
  PRINT *, 'dchi is zero'
  IF (outputdiags_flag) THEN
    PRINT *, 'Computing LS balance and energy diagnostics for t=0'
    CALL Calc_geost (LSfc(0))
    CALL Calc_hydro (LSfc(0), dims)
    CALL Energy (LSfc(0))
  END IF
  IF ((Vartype == 35) .OR. (Vartype == 4) .OR. forceNLmodel) THEN
    ! 3DFGAT or 4DVar, so compute the non-linear model trajectory
    PRINT *, '3DFGAT, 4DVar, or just propagating analysis throughout time window.'
    PRINT *, 'About to run nonlinear model for ', steps, ' time steps'
    CALL ABC_NL_ModelDriver_DA (LSfc(0:dasteps), dims, steps, dasteps)
    PRINT *, ' -- done'
    IF (outputdiags_flag) THEN
      DO t = 1, dasteps
        PRINT *, 'Computing LS balance diagnostics for t=', t
        CALL Calc_geost (LSfc(t))
        CALL Calc_hydro (LSfc(t), dims)
        CALL Energy (LSfc(t))
      END DO
    END IF
  ELSE
    ! 3DVar, so just copy t=0 state to all times
    PRINT *, '3DVar.  About to copy states'
    DO t = 1, dasteps
      LSfc(t) = LSfc(0)
    END DO
    PRINT *, ' -- done'
  END IF

  IF (outputdiags_flag) THEN
    ! Output the LS state throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/LS_Oloop', outerloop, '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting LS data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), LSfc(t), dims, dasteps+1, t, ModelStepsPerDAStep)
    END DO
  END IF


  ! Compute the reference model observations
  PRINT *, 'Computing model observations ...'
  CALL ModelObservations (dasteps, LSfc(0:dasteps), dims, dt_da, t0, &
                          Obs, .FALSE.)
  PRINT *, ' -- done'

  ! Set the remaining observation structure elements for this ref state and zero pert
  PRINT *, 'Setting model perts to zero ...'
  CALL ModelObservations_ZeroPert (Obs)
  PRINT *, ' -- done'


  IF (outputdiags_flag) THEN
    ! Output the observation structure
    WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/Obs_', outerloop, '_Iloop', innerloop, '.dat'
    PRINT *, 'Outputting obs structure ', TRIM(Diag_filename)
    CALL Write_Obs ( Diag_filename, 0, dt, 0, dt_da, 0, Obs )
  END IF

ELSE

  ! The reference values are already known
  PRINT *, 'dchi is non-zero'

  ! Convert the control variable pert to model space
  CALL U_trans (LSfc(0), dchi0, deltax(0), CVTdata, dims)

  IF (outputdiags_flag) THEN
    PRINT *, 'Computing pert balance diagnostics for t=0'
    CALL Calc_geost (deltax(0))
    CALL Calc_hydro (deltax(0), dims)
  END IF

  ! Propagate the initial perturbation throughout the time window
  IF ((VarType == 3) .OR. (VarType == 35)) THEN
    ! 3DVar or 3DFGAT, so initial perturbation does not change
    DO t = 1, dasteps
      deltax(t) = deltax(0)
    END DO
  ELSE
    ! 4DVar, so propagate perturbation
    ! Need to write tangent linear model
    ! Need to call balance diagnostics as above
    PRINT*, 'Error - 4DVar not yet implemented.'
    STOP
  END IF

!  IF (outputdiags_flag) THEN
!    ! Output the perturbation throughout the window
!    DO t = 0, dasteps
!      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/dx_Oloop', outerloop, &
!                            '_Iloop', innerloop, '.nc'
!      PRINT *, 'Time', t, ':  Outputting pert data to ', TRIM(Diag_filename)
!      CALL Write_state_2d (TRIM(Diag_filename), deltax(t), dims, dasteps+1, t, ModelStepsPerDAStep)
!    END DO
!  END IF

  ! Act with the linear observation operator
  CALL ModelObservations_linear (dasteps, LSfc(0:dasteps), deltax(0:dasteps), dims, &
                                 dt_da, t0, Obs )
    

END IF




! To do with the background term
! Compute the difference dchi0 - dchib0
CALL Initialise_CVs(diffcv, .FALSE.)
IF (.NOT.dchi0zero) THEN
  CALL Add_CVs(diffcv, dchi0)
END IF
IF (.NOT.dchib0zero) THEN
  CALL Subtract_CVs(diffcv, dchib0)
END IF



!IF (outputdiags_flag) THEN
!  ! Output the contribution to the gradient from the background term
!  WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/GradBg_Oloop', outerloop, &
!                        '_Iloop', innerloop, '.nc'
!  PRINT *, 'Outputting bg grad data to ', TRIM(Diag_filename)
!  CALL Write_state_2d (TRIM(Diag_filename), diffcv, dims, 1, 0, ModelStepsPerDAStep)
!END IF


IF (compute_grad) THEN

  ! ===================================
  ! --- ADJOINT PART OF CALCULATION ---
  ! ===================================

  ! Set the adjoint states to zero
  DO t = 0, dasteps
    CALL Initialise_model_vars (deltax_hat(t), .FALSE.)
  END DO

  ! Perform the adjoint of the model observations
  CALL ModelObservations_adj (dasteps, LSfc(0:dasteps), deltax_hat(0:dasteps), dims,  &
                              dt_da, t0, Obs, 1)
  ! At this stage the deltax_hat states are the Delta vectors in the documentation

  IF (outputdiags_flag .AND. (outerloop == 1) .AND. (innerloop == 1)) THEN
    ! Output the Delta throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/Delta_Oloop', outerloop, &
                            '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting Delta data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), deltax_hat(t), dims, dasteps+1, t, ModelStepsPerDAStep)
    END DO
  END IF

  ! Integrate the adjoint states backwards in time
  IF ((VarType == 3) .OR. (VarType == 35)) THEN
    ! 3DVar or 3DFGAT, so no time evolution of the adjoints
    DO t = dasteps-1, 0, -1
      CALL Add_model_vars(deltax_hat(t), deltax_hat(t+1), .TRUE.)
    END DO
  ELSE
    ! 4DVar, so propagate perturbation
    PRINT*, 'Error - 4DVar not yet implemented.'
    STOP
  END IF

  IF (outputdiags_flag .AND. (outerloop == 1) .AND. (innerloop == 1)) THEN
    ! Output the gradient of Jo throughout the window
    DO t = 0, dasteps
      WRITE (Diag_filename, '(A,A,I0.3,A,I0.3,A)') TRIM(outputdir), '/GradJo_Oloop', outerloop, &
                            '_Iloop', innerloop, '.nc'
      ! PRINT *, 'Time', t, ':  Outputting Delta data to ', TRIM(Diag_filename)
      CALL Write_state_2d (TRIM(Diag_filename), deltax_hat(t), dims, dasteps+1, t, ModelStepsPerDAStep)
    END DO
  END IF
  
  ! Operate with the adjoint of the CVT to obtain the gradient (observation contribution) in control space
  CALL U_trans_adj (LSfc(0), grad0, deltax_hat(0), CVTdata, dims)

  ! Add the background contribution to the gradient
  CALL Add_CVs(grad0, diffcv)


END IF


! ===================================
! --- VALUES OF COST FUNCTION ---
! ===================================

! The background part of the cost function
Jb_complex = InnerProdControlSpace (diffcv, diffcv, ComplexSpace=.TRUE., ignore_halos=.TRUE.) / 2.0
PRINT *, 'Jb = ', Jb_complex
Jb         = REAL(Jb_complex)

! The observation part of the cost function
Jo         = InnerProduct_obsself (Obs, .TRUE., 2)
PRINT *, 'Jo = ', Jo

! The total cost function
J          = Jb + Jo
PRINT *, 'J = ', J


END SUBROUTINE PenAndGrad
