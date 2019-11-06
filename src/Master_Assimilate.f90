PROGRAM Master_Assimilate

!*****************************************************
!*   Master code to do data assimilation with ABC    *
!*   model                                           *
!*                                                   *
!*   Ross Bannister, r.n.bannister@reading.ac.uk     *
!*   25/04/18                                        *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    Hybrid_opt,                  &
    Vartype,                     &
!    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    datadirAnal,                 &
    anal_file,                   &
    analinc_file,                &
    datadir_Bg,                  &
    Bg_file,                     &
    datadir_Obs,                 &
    Obs_file,                    &
    N_outerloops,                &
    N_innerloops_max,            &
    dt, dt_da, t0,               &
    mu, minus_mu,                &
    crit_inner,                  &
    dims_type,                   &
    ABC_type,                    &
    CV_type,                     &
    CVT_type,                    &
    diagnostics_file,            &
    Obs_type




IMPLICIT NONE

INCLUDE "InnerProdControlSpace.interface"
INCLUDE "Read_Obs.interface"
INCLUDE "DeAllocate_Obs.interface"
INCLUDE "PenAndGrad.interface"


! Declare variables
!==========================
CHARACTER(LEN=320)          :: Bg_filename, CVT_filename, Obs_filename
CHARACTER(LEN=320)          :: Scalar_diags_filename, Anal_filename, Analinc_filename
TYPE(ABC_type)              :: Bg0, dBg0, dx0, anal_inc
TYPE(ABC_type), ALLOCATABLE :: LSfc(:)
TYPE(dims_type)             :: dims
TYPE(CV_type)               :: dchiB0, dchi0_i, dchi0_ip1
TYPE(CV_type)               :: r_i, r_ip1, p_i, p_ip1, pert
TYPE(Obs_type), POINTER     :: Observations
INTEGER                     :: timesteps, DAtimesteps, maxtime
TYPE(CVT_type)              :: CVT
INTEGER                     :: k, i
LOGICAL                     :: dchiB0zero, dchi0zero, Converged_inner, ForceNLmodel
REAL(ZREAL8)                :: Jb, Jo, J, Jb_plus, Jo_plus, J_plus, Jb_minus, Jo_minus, J_minus
REAL(ZREAL8)                :: magnitude2_r_i, magnitude2_r_0, magnitude_adjust_inner, magnitude2_r_ip1
REAL(ZREAL8)                :: beta_i, alpha_i_i




PRINT*, '*************************************************************************'
PRINT*, 'Running Master_Assimilate'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions

  ! Main filenames (other files diagnostic files are used within the DA loops)
  ! Background state (input)
  Bg_filename             = TRIM(datadir_Bg)  // '/' // TRIM(Bg_file)
  ! CVT data (intput)
  CVT_filename            = TRIM(datadirCVT)  // '/' // TRIM(CVT_file)
  ! Obs data (input)
  Obs_filename            = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
  ! Diagnostics (output)
  Scalar_diags_filename   = TRIM(datadirAnal) // '/' // TRIM(diagnostics_file)
  ! Analysis (output)
  Anal_filename           = TRIM(datadirAnal) // '/' // TRIM(anal_file)
  ! Analysis increment (output)
  Analinc_filename        = TRIM(datadirAnal) // '/' // TRIM(analinc_file)

  OPEN (12, file=Scalar_diags_filename)
  WRITE (12,'(2A7, 8A15)') '#outer', 'inner', 'Jb', 'Jo', 'J', '|grad|', 'KE', 'BE', 'EE', 'TE'
  WRITE (12,'(2A7, 8A15)') '#-----', '-----', '--', '--', '-', '------', '--', '--', '--', '--'

  PRINT*, 'Reading in CVT data ...'
  CALL Read_Covs (CVT_filename, CVT,              &
                  .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
  PRINT *, '-- done'

  ! Read in the background state (-1 means the last item in the input file)
  PRINT*, 'Reading in background state ...'
  CALL Read_state_2d (Bg_filename, Bg0, dims, -1)
  PRINT*, '-- done'

  ! Set some commonly-used constants
  CALL Set_ht_dep_cons (dims)

  ! Compute balance and energy diagnostics of the background
  CALL Calc_geost (Bg0)
  CALL Calc_hydro (Bg0, dims)
  CALL Energy (Bg0)

  ! Initialise the observation structure (a linked list)
  NULLIFY(Observations)

  ! Read-in the observations
  PRINT*, 'Reading in observations ...'
  CALL Read_Obs (Observations, Obs_filename, dt, timesteps, dt_da, DAtimesteps, maxtime)
  PRINT *, 'Model timestep ', dt
  PRINT *, 'DA timeetep    ', dt_da
  PRINT *, '-- done'

  !Allocate the reference state array over the DA time steps
  ALLOCATE (LSfc(0:DAtimesteps))



  ! ===== Start the outer loop ================================================
  outerloop: DO k = 1, N_outerloops + 1

    PRINT *, '===== Starting outer loop', k


    ! Sort out the reference state and background perts
    PRINT *, 'Sorting out the initializations ...'
    IF (k == 1) THEN
      ! Set the reference state to the background
      LSfc(0) = Bg0
      ! This means that dBg0 and dchiB0 are both zero
      CALL Initialise_model_vars (dBg0, .FALSE.)
      CALL Initialise_CVs (dchiB0, .FALSE.)
      dchiB0zero = .TRUE.
    ELSE
      ! Compute the background as a perturbation
      dBg0 = Bg0
      CALL Subtract_model_vars (dBg0, LSfc(0), .TRUE.)
      ! Transform to control space
      CALL U_trans_inv (LSfc(0), dchiB0, dBg0, CVT, dims)
      dchiB0zero = .FALSE.
    END IF
    PRINT *, '-- done'

    ! In the inner loop, the first guess of the state is the reference
    ! This means that dchi0 is zero
    CALL Initialise_CVs (dchi0_i, .FALSE.)
    dchi0zero = .TRUE.

    ! Set to run the NL model in the PenAndGrad routine, despite assimilation mode
    ! This is used to run the NL model from the analysis (to find the next background) at the end of the DA
    ForceNLmodel = k > N_outerloops

    ! Compute the initial gradient of the cost function
    PRINT *, '===== Starting initial PenAndGrad routine'
    CALL PenAndGrad (VarType, DAtimesteps, timesteps, dims, t0,                    &
                     Observations, dchi0_i, dchi0zero, dchiB0, dchiB0zero,         &
                     CVT, LSfc, Jb, Jo, J, .TRUE., r_i,                            &
                     .TRUE., ForceNLmodel, datadirAnal, k, 0)
    PRINT *, '-- done'

    IF (k > N_outerloops) THEN
      WRITE (12,'(2I7, 8F15.5)') k, 0, Jb, Jo, J, 99999.0,                           &
                                 LSfc(0) % Kinetic_Energy, LSfc(0) % Buoyant_Energy, &
                                 LSfc(0) % Elastic_Energy, LSfc(0) % Total_Energy


    ELSE


      ! CALL Write_CV ('r_0_cv.nc', r_i, 14, CVT, dims % longs_v(1:nlongs), dims % full_levs(1:nlevs))


      ! The residual is minus this
      CALL Minus_CVs (r_i)
      magnitude2_r_i = REAL(InnerProdControlSpace(r_i, r_i, ComplexSpace=.TRUE.))
      ! Store this value as the initial magnitude squared of the residual
      magnitude2_r_0 = magnitude2_r_i

      ! Output initial diagnostics
      WRITE (12,'(2I7, 8F15.5)') k, 0, Jb, Jo, J, SQRT(magnitude2_r_i),              &
                                 LSfc(0) % Kinetic_Energy, LSfc(0) % Buoyant_Energy, &
                                 LSfc(0) % Elastic_Energy, LSfc(0) % Total_Energy


      ! The initial search direction is the initial residual
      p_i = r_i

      ! ===== Start the inner loop ================================================
      i = 1
      innerloop: DO

        PRINT *, 'Starting inner loop:', i

        ! Do a line minimization in the search direction
        ! Compute J at a positively perturbed state
        CALL Add_pert_CVs (pert, dchi0_i, p_i, mu)
        PRINT *, 'Starting PenAndGrad routine for J+'
        CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0,            &
                          Observations, pert, .FALSE., dchiB0, dchiB0zero,      &
                          CVT, LSfc, Jb_plus, Jo_plus, J_plus, .FALSE., r_i,    &
                          .FALSE., .FALSE., datadirAnal, k, i )
        PRINT *, 'Positive pert (Jb+, Jo+, J+):', Jb_plus, Jo_plus, J_plus


        ! Compute J at a negatively perturbed state
        CALL Add_pert_CVs (pert, dchi0_i, p_i, minus_mu)
        PRINT *, 'Starting PenAndGrad routine for J-'
        CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0,            &
                          Observations, pert, .FALSE., dchiB0, dchiB0zero,      &
                          CVT, LSfc, Jb_minus, Jo_minus, J_minus, .FALSE., r_i, &
                          .FALSE., .FALSE., datadirAnal, k, i )
        PRINT *, 'Negative pert (Jb-, Jo-, J-):', Jb_minus, Jo_minus, J_minus

        beta_i = mu * (J_minus - J_plus) / (2.0 * (J_minus + J_plus - 2.0 * J))

        ! The adjustment to the state is beta_i * p_i
        ! What is the size of this adjustment?
        magnitude_adjust_inner = beta_i * SQRT(REAL(InnerProdControlSpace(p_i, p_i, ComplexSpace=.TRUE.)))
        !PRINT *, 'Change in dchi:', k, i, magnitude_adjust_inner

        ! Update the state
        CALL Add_pert_CVs (dchi0_ip1, dchi0_i, p_i, beta_i)
        dchi0zero = .FALSE.

        ! Compute a new gradient
        PRINT *, 'Starting main PenAndGrad routine for inner loop', i
        CALL PenAndGrad ( VarType, DAtimesteps, timesteps, dims, t0,              &
                          Observations, dchi0_ip1, dchi0zero, dchiB0, dchiB0zero, &
                          CVT, LSfc, Jb, Jo, J, .TRUE., r_ip1,                    &
                          .TRUE., .FALSE., datadirAnal, k, i )

        ! The residual is minus this
        CALL Minus_CVs (r_ip1)
        magnitude2_r_ip1 = REAL(InnerProdControlSpace(r_ip1, r_ip1, ComplexSpace=.TRUE.))

        ! Output diagnostics
        WRITE (12,'(2I7, 8F15.5)') k, i, Jb, Jo, J, SQRT(magnitude2_r_i),         &
                                   0.0, 0.0, 0.0, 0.0


        !PRINT *, 'Next residual:', k, i, SQRT(magnitude2_r_ip1)

        ! Have we converged?
        Converged_inner = ( SQRT(magnitude2_r_ip1 / magnitude2_r_0) < crit_inner )

        IF (.NOT.Converged_inner) THEN
          alpha_i_i = magnitude2_r_ip1 / magnitude2_r_i

          ! Compute the new search direction
          CALL Add_pert_CVs (p_ip1, r_ip1, p_i, alpha_i_i)

          ! Increment inner loop counter, and shift variables
          i              = i + 1
          r_i            = r_ip1
          magnitude2_r_i = magnitude2_r_ip1
          p_i            = p_ip1
          dchi0_i        = dchi0_ip1

          IF (i > N_innerloops_max) THEN
            EXIT innerloop
          END IF
        ELSE
          EXIT innerloop
        END IF
      END DO innerloop

      ! Compute the increment in model space
      CALL U_trans (LSfc(0), dchi0_ip1, dx0, CVT, dims)

      ! Add to reference state
      CALL Add_model_vars (LSfc(0), dx0, .TRUE.)

    END IF

  END DO outerloop


  ! Compute the overall analysis increment
  anal_inc = LSfc(0)
  CALL Subtract_model_vars (anal_inc, Bg0, .FALSE.)

  ! Output the analysis and analysis increment
  PRINT*, 'Writing the analysis ...'
  CALL Write_state_2d (Anal_filename, LSfc(0), dims, 1, 0, 0)
  PRINT*, '  -- done'
  PRINT*, 'Writing the analysis increment ...'
  CALL Write_state_2d (Analinc_filename, anal_inc, dims, 1, 0, 0)
  PRINT*, '  -- done'

  CLOSE (12)

  ! Deallocate
  PRINT*, 'Deallocating ...'
  CALL DeAllocate_Obs (Observations)
  DEALLOCATE (LSfc)
  PRINT*, '  -- done'

END PROGRAM Master_Assimilate
