PROGRAM Master_Calibration

!*****************************************************
!*   Code to do either of the following              *
!*                                                   *
!*   Stage 1 UM to ABC forecasts                     *
!*                                                   *
!*   Stage 2 compute ABC perturbations               *
!*                                                   *
!*   Stage 3 Compute balanced mass regression        *
!*                                                   *
!*   Stage 4 param transform (convert model perts to *
!*                        parameters)                *
!*   Stage 5 calibration (calibrate spatial          *
!*                        transforms)                *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    zero,                        &
    nlongs, nlevs,               &
    dims_type,                   &
    UM_type,                     &
    ABC_type,                    &
    CV_type,                     &
    CVT_type,                    &
    datadirCVT,                  &
    CVT_file,                    &
    diagnostics_file,            &
    dt, deltat,                  &
    ntimesteps,                  &
    NEnsmax, NEns,               &
    EnsDirs,                     &
    NEnsMems,                    &
    NNMCmax, NNMC,               &
    NMCDirs,                     &
    Nlatsmax, Nlats,             &
    latindex,                    &
    datadirABCfcs,               &
    datadirABCperts,             &
    datadirRegression,           &
    datadirConParams,            &
    CalibRunStage,               &
    Regular_vert_grid,           &
    Adv_tracer,                  &
    CVT,                         &
    VertSmoothPoints,            &
    HorizSmoothPoints,           &
    ForceCor,                    &
    sr_nlongs


IMPLICIT NONE

INCLUDE "Boundaries.interface"
INCLUDE "Boundaries_CV.interface"
INCLUDE "Write_Covs.interface"


! Declare variables
!==========================
INTEGER                     :: ens, mem, lat, lev, long, item, t, Neffmems, Nloops, Nstates, loop, k
REAL(ZREAL8)                :: Neffmems_r, NTotalStates_r, population_r, Nens_r, alpha
CHARACTER(LEN=320)          :: state_file
CHARACTER(LEN=320)          :: paramfile
CHARACTER(LEN=320)          :: meanfile
CHARACTER(LEN=320)          :: ABCfile
CHARACTER(LEN=320)          :: CVT_filename
TYPE(UM_type)               :: um_data
TYPE(dims_type)             :: dims
TYPE(ABC_type)              :: ABC_data_a, ABC_data_b
TYPE(ABC_type)              :: ABC_mean
REAL(ZREAL8), ALLOCATABLE   :: r_total(:,:,:)
REAL(ZREAL8), ALLOCATABLE   :: r_bal(:,:,:)
REAL(ZREAL8), ALLOCATABLE   :: Inv_Cov_rbalrbal(:,:)
REAL(ZREAL8), ALLOCATABLE   :: Cov_temp(:,:)
TYPE(CV_type)               :: ControlVar, ControlVar_mean, ControlVar_1


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_Calibration'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions


  IF ((Nens <= 0) .AND. (NNMC <= 0)) THEN
    PRINT *, 'One or more of Nens and NNMC must be positive'
    STOP
  END IF

  IF ((Nlats < 1) .OR. (Nlats > Nlatsmax)) THEN
    PRINT *, 'Nlats must be between 1 and ', Nlatsmax
    PRINT *, 'Either run with different Nlats, or recompile changing Nlatsmax'
    STOP
  END IF

  IF (Nens > 0) THEN
    ! We want to calibrate from ensembles
    PRINT *, 'Calibrating with ensembles'
    IF (Nens > Nensmax) THEN
      PRINT *, 'Error: Nens needs to be smaller than or equal to Nensmax'
      PRINT *, 'You need to recompile with a larger value of Nensmax'
      STOP
    END IF
  END IF
  IF (NNMC > 0) THEN
    ! We want to calibrate from NMC data
    PRINT *, 'Calibrating with NMC pairs'
    IF (NNMC > NNMCmax) THEN
      PRINT *, 'Error: NNMC needs to be smaller than or equal to NNMCmax'
      PRINT *, 'You need to recompile with a larger value of NNMCmax'
      STOP
    END IF
  END IF

  IF ( (dt/deltat /= 2.0 )) THEN
    PRINT*, 'Error: dt /deltat .NE. 2'
    STOP
  ENDIF





  SELECT CASE (CalibRunStage)

  ! ====================================================================
  CASE (1)
    PRINT *, 'Stage 1 of the calibration - generate ABC-style forecasts'


    IF (Nens > 0) THEN
      ! We want to calibrate the CVT from ensembles
      CALL Initialise_um_data (um_data)
      CALL Initialise_dims (dims)
      CALL Initialise_model_vars (ABC_data_a, .FALSE.)
      CALL Initialise_model_vars (ABC_data_b, .FALSE.)

      DO ens = 1, Nens
        item = 0
        DO mem = 1, NEnsMems
          WRITE (state_file, '(A,A,I0.3,A)') TRIM(EnsDirs(ens)), '/Member', mem, '.nc'
          PRINT *, TRIM(state_file)
          DO lat = 1, Nlats
            item = item + 1
            PRINT *, 'Lat index = ', latindex(lat)

            ! Read in raw UM data, store in init_um_data
            PRINT *, 'Reading UM data ...'
            CALL Read_um_data_2d (um_data, state_file, latindex(lat))
            PRINT *, '--done'

            ! Process um data
            ! Set grid
            IF ((ens == 1) .AND. (mem == 1) .AND. (lat == 1)) THEN
              PRINT *, 'Setting Grid ...'
              CALL Set_grid (um_data, dims, Regular_vert_grid)
              CALL Set_ht_dep_cons (dims)
              PRINT *, '--done'
            END IF

            ! Define variables for simplified model from UM data, store in init_state
            PRINT *, 'Processing UM data to make it compatible with simplified model ...'
            CALL Process_um_data (um_data, ABC_data_a, dims)
            PRINT *, '--done'

            ! Apply boundary conditions
            CALL Boundaries (ABC_data_a, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_r=.TRUE., &
                             set_b=.TRUE., set_rho=.TRUE., set_tracer=Adv_tracer)

            ! Output initial conditions (to check working)
            !WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCfcs), '/IC_Ens', ens, '_Item', item, '.nc'
            !CALL Write_state_2d (testfile, ABC_data_a, dims, 1, 0, 0)

            ! Integrate this model forward in time
            DO t = 1, ntimesteps
              CALL ABC_NL_model (ABC_data_a, ABC_data_b, dims)
              ABC_data_a = ABC_data_b
            END DO

            ! Output short forecast
            WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCfcs), '/FC_Ens', ens, '_Item', item, '.nc'
            CALL Write_state_2d (ABCfile, ABC_data_a, dims, 1, 0, 0)

          END DO
        END DO
      END DO


    ELSE
      ! We want to calibrate the CVT from NMC data
      PRINT *, 'NMC method not yet implemented'

    END IF


    ! Create covariance file, and put options and A, B, C, f into it
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    PRINT *, 'Creating CVT file: ', TRIM(CVT_filename)
    CALL Write_Covs (CVT_filename, CVT,                     &
                     dims % longs_v(1:nlongs),              &
                     dims % half_levs(1:nlevs),             &
                     newfile = .TRUE., output_params = .TRUE.)



  ! ====================================================================
  CASE (2)
    PRINT *, 'Stage 2 of the calibration - compute perturbations'
    IF (Nens > 0) THEN
      CALL Initialise_dims (dims)
      CALL Initialise_model_vars (ABC_data_a, .FALSE.)
      Neffmems   = NEnsMems * Nlats
      Neffmems_r = REAL(Neffmems)

      PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
      DO ens = 1, Nens

        ! Compute the mean for this ensemble
        CALL Initialise_model_vars (ABC_mean, .FALSE.)
        item = 0
        DO mem = 1, NEnsMems
          DO lat = 1, Nlats
            item = item + 1
            ! Read-in a state found in stage 1
            WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCfcs), '/FC_Ens', ens, '_Item', item, '.nc'
            PRINT *, 'Reading file ', TRIM(ABCfile)
            CALL Read_state_2d (ABCfile, ABC_data_a, dims, 1)
            PRINT *, '-- done'
            IF ((ens == 1) .AND. (mem == 1) .AND. (lat == 1)) THEN
              CALL Set_ht_dep_cons (dims)
            END IF

            CALL Add_model_vars(ABC_mean, ABC_data_a, .TRUE.)

          END DO
        END DO

        CALL Div_model_cons(ABC_mean, Neffmems_r, .TRUE.)

        ! Write-out the ensemble mean
        WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
        PRINT *, 'Writing mean state ', TRIM (meanfile)
        CALL Write_state_2d (meanfile, ABC_mean, dims, 1, 0, 0)
        PRINT *, '-- done'


        ! Re-read forecasts, and compute and then output the perts
        item = 0
        DO mem = 1, NEnsMems
          DO lat = 1, Nlats
            item = item + 1
            ! Read-in a state found in stage 1
            WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCfcs), '/FC_Ens', ens, '_Item', item, '.nc'
            PRINT *, 'Reading file ', TRIM(ABCfile)
            CALL Read_state_2d (ABCfile, ABC_data_a, dims, 1)
            PRINT *, '-- done'

            CALL Subtract_model_vars(ABC_data_a, ABC_mean, .TRUE.)

            WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
            PRINT *, 'Writing perturbation file ', TRIM(ABCfile)
            CALL Write_state_2d (ABCfile, ABC_data_a, dims, 1, 0, 0)
            PRINT *, '-- done'

          END DO
        END DO

      END DO


    ELSE
      ! We want to calibrate the CVT from NMC data
      PRINT *, 'NMC method not yet implemented'

    END IF








  ! ====================================================================
  CASE (3)
    PRINT *, 'Stage 3 of the calibration - compute regression matrix (if required)'
    PRINT *, 'Reading options from CVT file'
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    CALL Read_Covs (CVT_filename, CVT,           &
                    .TRUE.,                      &
                    .FALSE., .FALSE.,            &
                    .FALSE., .FALSE.)
    PRINT *, '-- done'

    ! Initialize some arrays
    CVT % Cov_rtotrbal(1:nlevs,1:nlevs) = 0.0
    CVT % Cov_rbalrbal(1:nlevs,1:nlevs) = 0.0
    CVT % Regression(1:nlevs,1:nlevs)   = 0.0


    IF ((CVT % CVT_param_opt_reg == 1) .AND. (CVT % CVT_param_opt_gb < 3)) THEN
      IF (Nens > 0) THEN

        CALL Initialise_dims (dims)
        Neffmems   = NEnsMems * Nlats
        Neffmems_r = REAL(Neffmems)
        Nens_r     = REAL(Nens)
        ALLOCATE (r_total(1:nlongs,1:nlevs, 1:Neffmems))
        ALLOCATE (r_bal(1:nlongs,1:nlevs, 1:Neffmems))
        ALLOCATE (Inv_Cov_rbalrbal(1:nlevs,1:nlevs))
        ALLOCATE (Cov_temp(1:nlevs,1:nlevs))

      	DO ens = 1, Nens
          item = 0

          DO mem = 1, NEnsMems
            DO lat = 1, Nlats
              item = item + 1

              ! Read-in a perturbation state found in stage 2
              WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
              PRINT *, 'Reading file ', TRIM(ABCfile)
              CALL Read_state_2d (ABCfile, ABC_data_a, dims, 1)
              PRINT *, '-- done'

              ! Extract the total r field and store for computation
              r_total(1:nlongs,1:nlevs, item) = ABC_data_a % r(1:nlongs,1:nlevs)

              IF ((ens == 1) .AND. (item == 1)) THEN
                ! Set dimensions
                CALL Set_ht_dep_cons (dims)
                ! Output r field as a sample
                state_file = TRIM(datadirRegression) // '/r_001.nc'
                CALL Write_one_field (state_file, nlongs, nlevs, &
                                      r_total(1:nlongs,1:nlevs,1), 'r_total')
              END IF

              ! Compute the balanced r field and store for computation
              ! First compute psi and chi from u and v (psi is the important one)
              PRINT *, 'Computing psi and chi'
              CALL Helmholtz_inv (ControlVar % v1(1:nlongs,1:nlevs),     &
                		  ControlVar % v2(1:nlongs,1:nlevs),     &
                		  ABC_data_a % u(0:nlongs+1,1:nlevs),    &
                		  ABC_data_a % v(0:nlongs+1,1:nlevs))

              IF ((ens ==1) .AND. (item ==1)) THEN
                ! Output streamfunction as a sample
                state_file = TRIM(datadirRegression) // '/psi_001.nc'
                CALL Write_one_field (state_file, nlongs, nlevs, &
                                      ControlVar % v1(1:nlongs,1:nlevs), 'psi')
              END IF

              PRINT *, '-- done'
              CALL Boundaries_CV (ControlVar, set_1=.TRUE., set_2=.TRUE.)

              IF (CVT % CVT_param_opt_gb == 1) THEN
                ! Analytical balance
                CALL LinearBal_r (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                                  r_bal(1:nlongs, 1:nlevs, item))

                IF ((ens ==1) .AND. (item ==1)) THEN
                  ! Output balanced r as a sample
                  state_file = TRIM(datadirRegression) // '/rbal_001.nc'
                  CALL Write_one_field (state_file, nlongs, nlevs, &
                                        r_bal(1:nlongs,1:nlevs,1), 'r_bal')
                END IF

      	      ELSE
      		      ! Statistical balance
                PRINT *, 'Statistical mass/wind balance not implemented'
                STOP
      	      END IF
            END DO
          END DO

          ! Compute vertical covariances of total and balanced pressures
          PRINT *, 'Computing vertical covariance matrices'
          PRINT *, 'There are ', Neffmems, ' effective members'
          CALL Calc_vert_cov1 (Neffmems,                                 &
                               r_bal(1:nlongs,1:nlevs,1:Neffmems),       &
                               r_bal(1:nlongs,1:nlevs,1:Neffmems),       &
                               Cov_temp(1:nlevs,1:nlevs))
          PRINT *, '-- done 1'
          CVT % Cov_rbalrbal(1:nlevs,1:nlevs) = CVT % Cov_rbalrbal(1:nlevs,1:nlevs) +&
                                                Cov_temp(1:nlevs,1:nlevs)
          CALL Calc_vert_cov1 (Neffmems,                                 &
                               r_total(1:nlongs,1:nlevs,1:Neffmems),     &
                               r_bal(1:nlongs,1:nlevs,1:Neffmems),       &
                               Cov_temp(1:nlevs,1:nlevs))
          PRINT *, '-- done 2'
          CVT % Cov_rtotrbal(1:nlevs,1:nlevs) = CVT % Cov_rtotrbal(1:nlevs,1:nlevs) +&
                                                Cov_temp(1:nlevs,1:nlevs)
        END DO

        ! Normalize the covariance matrices
        CVT % Cov_rbalrbal(1:nlevs,1:nlevs) = CVT % Cov_rbalrbal(1:nlevs,1:nlevs) / Nens_r
        CVT % Cov_rtotrbal(1:nlevs,1:nlevs) = CVT % Cov_rtotrbal(1:nlevs,1:nlevs) / Nens_r

        ! Find the inverse covariance of the balanced covariances
        Inv_Cov_rbalrbal(1:nlevs,1:nlevs) = CVT % Cov_rbalrbal(1:nlevs,1:nlevs)
        CALL InverseSymMat (nlevs, Inv_Cov_rbalrbal)

        ! Form the vertical mass regression matrix
        CVT % Regression(1:nlevs,1:nlevs) = MATMUL(CVT % Cov_rtotrbal, Inv_Cov_rbalrbal)

        DEALLOCATE (r_total, r_bal, Inv_Cov_rbalrbal, Cov_temp)

      ELSE
        ! We want to calibrate the CVT from NMC data
        PRINT *, 'NMC method not yet implemented'
        STOP
      END IF
    ELSE
      PRINT *, 'There is no requirement/request to use/compute a regression matrix'
    END IF

    ! Output this regression matrix for use later in stage 4
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    PRINT *, 'Outputting regression data in CVT file: ', TRIM(CVT_filename)
    CALL Write_Covs (CVT_filename, CVT,                     &
                     dims % longs_v(1:nlongs),              &
                     dims % half_levs(1:nlevs),             &
                     output_regression = .TRUE.)





  ! ====================================================================
  CASE (4)
    PRINT *, 'Stage 4 of the calibration - perform parameter transform on perts'
    ! Read-in parameters, options, and vertical regression matrix
    PRINT *, 'Reading regression data'
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    CALL Read_Covs (CVT_filename, CVT,           &
                    .TRUE.,                      &
                    .FALSE., .FALSE.,            &
                    .FALSE., .TRUE.)
    PRINT *, '-- done'


    PRINT *, 'Cov options:'
    PRINT *, '------------'
    PRINT *, CVT % CVT_order, CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, &
             CVT % CVT_param_opt_ab, CVT % CVT_param_opt_reg
    PRINT *, '------------'


    IF (Nens > 0) THEN
      CALL Initialise_dims (dims)
      Neffmems   = NEnsMems * Nlats
      Neffmems_r = REAL(Neffmems)

      PRINT *, 'There are ', Neffmems, ' effective ensemble members.'
      DO ens = 1, Nens

        ! Read-in the ensemble mean for this ensemble
        WRITE (meanfile, '(A,A,I0.3,A)') TRIM(datadirABCperts), '/MeanABC', ens, '.nc'
        PRINT *, 'Reading mean state ', TRIM (meanfile)
        CALL Read_state_2d (meanfile, ABC_mean, dims, 1)
        PRINT *, '-- done'

        ! Read-in the ABC perts and convert to parameter perts
        item = 0
        DO mem = 1, NEnsMems
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in a pert found from stage 2
            WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
            PRINT *, 'Reading file ', TRIM(ABCfile)
            CALL Read_state_2d (ABCfile, ABC_data_a, dims, 1)
            PRINT *, '-- done'
            IF ((ens == 1) .AND. (item == 1)) THEN
              CALL Set_ht_dep_cons (dims)
            END IF

            ! Pass through the inverse parameter transform
            CALL U_p_inv (ABC_mean, ControlVar, ABC_data_a, CVT % CVT_order,                        &
                          CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,   &
                          CVT % CVT_param_opt_reg,                                                  &
                          CVT % Regression(1:nlevs,1:nlevs), dims,                                  &
                          ((ens == 1) .AND. (item ==1)),  &  ! Output diagnostics only if this is met
                          datadirConParams)                  ! For diagnostic location

            ! Write parameter pert
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', ens, '_Item', item, '.nc'
            PRINT *, 'Writing parameters to file ', TRIM (paramfile)
            CALL Write_CV (paramfile, ControlVar, 1, CVT, dims % longs_v(1:nlongs), dims % half_levs(1:nlevs))


          END DO
        END DO

      END DO


    ELSE
      ! We want to calibrate the CVT from NMC data
      PRINT *, 'NMC method not yet implemented'
      STOP
    END IF





  ! ====================================================================
  CASE (5)
    PRINT *, 'Stage 5 of the calibration - calibrate stddevs and spatial transforms'

    PRINT *, 'Reading options from CVT file'
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    CALL Read_Covs (CVT_filename, CVT,           &
                    .TRUE.,                      &
                    .FALSE., .FALSE.,            &
                    .FALSE., .FALSE.)
    PRINT *, '-- done'

    ! Read-in a sample model state to get dimension information
    WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', 1, '_Item', 1, '.nc'
    PRINT *, 'Reading file (for dimension information) ', TRIM(ABCfile)
    CALL Read_state_2d (ABCfile, ABC_data_a, dims, 1)
    PRINT *, '-- done'
    CALL Set_ht_dep_cons (dims)


    IF (Nens > 0) THEN
      ! Perts have been generated with the ensemble method
      ! Call the number of ensembles now Nloops
      Nloops   = Nens
      Nstates  = NEnsMems
      Neffmems = NEnsMems * Nlats
    ELSE
      ! Perts have been generated with the NMC method
      Nloops   = NNMC
      Nstates  = 1
      Neffmems = Nlats
    END IF

    Neffmems_r     = REAL(Neffmems)
    NTotalStates_r = REAL(Neffmems * Nloops)




    ! --- Compute the standard deviations of the parameters ---

    PRINT *, 'Computing the standard deviations'
    ! Initialise the master running total

    DO loop = 1, Nloops      ! For ensembles, this loop is over ensembles
      item = 0
      ! Initialise the running total
      CALL Initialise_CVs (ControlVar_mean, .FALSE.)
      DO mem = 1, Nstates    ! For ensembles, this loops is over ensemble members
        DO lat = 1, Nlats
          item = item + 1

          ! Read-in this perturbation
          WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
          CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

          ! Increment the running total
          CALL Add_CVs (ControlVar_mean, ControlVar)

        END DO
      END DO

      ! Find the mean for this set of perturbations
      CALL Div_CV_cons (ControlVar_mean, Neffmems_r)

      ! Repeat the reading-in to find the standard deviations
      item = 0
      DO mem = 1, Nstates
        DO lat = 1, Nlats
          item = item + 1

          ! Read-in this perturbation
          WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
          CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

          ! Remove the mean found above
          CALL Subtract_CVs (ControlVar, ControlVar_mean)

          ! Compute contribution to the standard deviations
          CALL CVT_Calibration_stddevs (CVT, ControlVar, zero)  !0.0 means contribute to stddevs

        END DO
      END DO

    END DO

    ! Square-root and normalize for standard deviations, then do required averaging as CVT_stddev_opt
    CALL CVT_Calibration_stddevs (CVT, ControlVar, NTotalStates_r)  !ControlVar is not used in this call

    ! Do any smoothing that is requested
    IF (VertSmoothPoints > 0) THEN
      ! Vertical smoothing of standard deviations
      PRINT *, 'Doing smoothing of the standard deviations in the vertical'
      DO long = 1, nlongs
        CALL Smooth (nlevs, CVT % sigma1(long, 1:nlevs), VertSmoothPoints)
        CALL Smooth (nlevs, CVT % sigma2(long, 1:nlevs), VertSmoothPoints)
        CALL Smooth (nlevs, CVT % sigma3(long, 1:nlevs), VertSmoothPoints)
        CALL Smooth (nlevs, CVT % sigma4(long, 1:nlevs), VertSmoothPoints)
        CALL Smooth (nlevs, CVT % sigma5(long, 1:nlevs), VertSmoothPoints)
        CALL Smooth (nlevs, CVT % sigma6(long, 1:nlevs), VertSmoothPoints)
      END DO
    END IF
    IF (HorizSmoothPoints > 0) THEN
      PRINT *, 'Doing smoothing of the standard deviations in the horizontal'
      DO lev = 1, nlevs
        CALL Smooth (nlongs, CVT % sigma1(1:nlongs, lev), HorizSmoothPoints)
        CALL Smooth (nlongs, CVT % sigma2(1:nlongs, lev), HorizSmoothPoints)
        CALL Smooth (nlongs, CVT % sigma3(1:nlongs, lev), HorizSmoothPoints)
        CALL Smooth (nlongs, CVT % sigma4(1:nlongs, lev), HorizSmoothPoints)
        CALL Smooth (nlongs, CVT % sigma5(1:nlongs, lev), HorizSmoothPoints)
        CALL Smooth (nlongs, CVT % sigma6(1:nlongs, lev), HorizSmoothPoints)
      END DO
    END IF


    ! Write-out the standard deviations
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    PRINT *, 'Outputting stddevs in CVT file: ', TRIM(CVT_filename)
    CALL Write_Covs (CVT_filename, CVT,                     &
                     dims % longs_v(1:nlongs),              &
                     dims % half_levs(1:nlevs),             &
                     output_stddevs = .TRUE.)



    ! --- Compute the spatial transforms (this depends on the order of the transforms) ---

    SELECT CASE (CVT % CVT_order)


    CASE (1)
      ! The original MetO order of transforms
      ! -------------------------------------

      PRINT *, 'Computing the vertical transform (original MetO transform order)'

      DO loop = 1, Nloops
        PRINT *, '  Finding the mean, loop', loop
        item = 0
        ! Initialise the running total
        CALL Initialise_CVs (ControlVar_mean, .FALSE.)
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Increment the running total
            CALL Add_CVs (ControlVar_mean, ControlVar)
          END DO
        END DO

        ! Find the mean for this set of perturbations
        CALL Div_CV_cons (ControlVar_mean, Neffmems_r)

        PRINT *, '  Finding the vert covs loop', loop
        ! Repeat the reading-in to find the standard deviations
        item = 0
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Remove the mean found above
            CALL Subtract_CVs (ControlVar, ControlVar_mean)

            ! Compute contribution to the vertical covs
            CALL CVT_Calibration_vertcovs (CVT, ControlVar, zero)  !0.0 means contribute to vert covs

          END DO
        END DO

      END DO

      population_r = REAL(nlongs) * NTotalStates_r
      ! Normalize for vertical covariances
      CALL CVT_Calibration_vertcovs (CVT, ControlVar, population_r)  !ControlVar is not used in this call

      ! Due to approximations made, the vertical matrix will not precisely be a correlation matrix
      ! Scale to make sure that it is a correlation matrix
      IF (ForceCor) THEN
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode1(1:nlevs,1:nlevs,1))
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode2(1:nlevs,1:nlevs,1))
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode3(1:nlevs,1:nlevs,1))
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode4(1:nlevs,1:nlevs,1))
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode5(1:nlevs,1:nlevs,1))
        CALL Ensure_correlation_matrix (nlevs, CVT % VertMode6(1:nlevs,1:nlevs,1))
      END IF

      ! Now compute the vertical eigenvectors and eigenvalues
      CALL VertEigens (CVT)

      ! Output
      CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
      PRINT *, 'Outputting vertical covariances in CVT file: ', TRIM(CVT_filename)
      CALL Write_Covs (CVT_filename, CVT,                     &
                       dims % longs_v(1:nlongs),              &
                       dims % half_levs(1:nlevs),             &
                       output_vert = .TRUE.)




      PRINT *, 'Computing the horizontal transform (original MetO transform order)'

      DO loop = 1, Nloops
        PRINT *, '  Finding the mean, loop', loop
        item = 0
        ! Initialise the running total
        CALL Initialise_CVs (ControlVar_mean, .FALSE.)
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Perform the inverse vertical transform
            CALL U_v_inv(ControlVar_1, ControlVar, CVT)

            ! Do FFT
            CALL U_h_inv(ControlVar, ControlVar_1, CVT, .TRUE.)

            ! Increment the running total
            CALL Add_CVs (ControlVar_mean, ControlVar)
          END DO
        END DO

        ! Find the mean for this set of perturbations
        CALL Div_CV_cons (ControlVar_mean, Neffmems_r)

        PRINT *, '  Finding the horiz spectra, loop', loop
        ! Repeat the reading-in to find the standard deviations
        item = 0
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Perform the inverse vertical transform
            CALL U_v_inv(ControlVar_1, ControlVar, CVT)

            ! Do FFT
            CALL U_h_inv(ControlVar, ControlVar_1, CVT, .TRUE.)

            ! Remove the mean found above
            CALL Subtract_CVs (ControlVar, ControlVar_mean)

            ! Compute contribution to the horizontal variances
            CALL CVT_Calibration_horizcovs (CVT, ControlVar, zero)     !0.0 means contribute to horiz covs

          END DO
        END DO

      END DO

      ! Normalize and square-root for standard deviations of horizontal spectra
      CALL CVT_Calibration_horizcovs (CVT, ControlVar, NTotalStates_r)  !ControlVar is not used in this call

      IF (ForceCor) THEN
        ! Modify to ensure that the horizontal transform represents a correlation matrix
        DO lev = 1, nlevs
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV1(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV1(1:nlongs/2+1,lev)))
          CVT % HorizEV1(1:nlongs/2+1,lev) = alpha * CVT % HorizEV1(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV2(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV2(1:nlongs/2+1,lev)))
          CVT % HorizEV2(1:nlongs/2+1,lev) = alpha * CVT % HorizEV2(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV3(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV3(1:nlongs/2+1,lev)))
          CVT % HorizEV3(1:nlongs/2+1,lev) = alpha * CVT % HorizEV3(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV4(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV4(1:nlongs/2+1,lev)))
          CVT % HorizEV4(1:nlongs/2+1,lev) = alpha * CVT % HorizEV4(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV5(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV5(1:nlongs/2+1,lev)))
          CVT % HorizEV5(1:nlongs/2+1,lev) = alpha * CVT % HorizEV5(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV6(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV6(1:nlongs/2+1,lev)))
          CVT % HorizEV6(1:nlongs/2+1,lev) = alpha * CVT % HorizEV6(1:nlongs/2+1,lev)
        END DO
      END IF

      ! Output
      CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
      PRINT *, 'Outputting horizontal spectra in CVT file: ', TRIM(CVT_filename)
      CALL Write_Covs (CVT_filename, CVT,                     &
                       dims % longs_v(1:nlongs),              &
                       dims % half_levs(1:nlevs),             &
                       output_horiz = .TRUE.)







    CASE (2)
      ! The reversed MetO order of transforms
      ! -------------------------------------

      PRINT *, 'Computing the horizontal transform (reversed MetO transform order)'

      DO loop = 1, Nloops
        PRINT *, '  Finding the mean, loop', loop
        item = 0
        ! Initialise the running total
        CALL Initialise_CVs (ControlVar_mean, .FALSE.)
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Do FFT
            CALL U_h_inv(ControlVar_1, ControlVar, CVT, .TRUE.)

            ! Increment the running total
            CALL Add_CVs (ControlVar_mean, ControlVar_1)

          END DO
        END DO

        ! Find the mean for this set of perturbations
        CALL Div_CV_cons (ControlVar_mean, Neffmems_r)

        ! Repeat the reading-in to find the standard deviations
        PRINT *, '  Finding the horizontal spectra, loop', loop
        item = 0
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Do FFT
            CALL U_h_inv(ControlVar_1, ControlVar, CVT, .TRUE.)

            ! Remove the mean found above
            CALL Subtract_CVs (ControlVar_1, ControlVar_mean)

            ! Compute contribution to the horizontal variances
            CALL CVT_Calibration_horizcovs (CVT, ControlVar_1, zero)     !0.0 means contribute to horiz covs

          END DO
        END DO

      END DO

      ! Normalize and square-root for standard deviations of horizontal spectra
      CALL CVT_Calibration_horizcovs (CVT, ControlVar, NTotalStates_r)  !ControlVar is not used in this call

      IF (ForceCor) THEN
        ! Modify to ensure that the horizontal transform represents a correlation matrix
        DO lev = 1, nlevs
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV1(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV1(1:nlongs/2+1,lev)))
          CVT % HorizEV1(1:nlongs/2+1,lev) = alpha * CVT % HorizEV1(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV2(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV2(1:nlongs/2+1,lev)))
          CVT % HorizEV2(1:nlongs/2+1,lev) = alpha * CVT % HorizEV2(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV3(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV3(1:nlongs/2+1,lev)))
          CVT % HorizEV3(1:nlongs/2+1,lev) = alpha * CVT % HorizEV3(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV4(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV4(1:nlongs/2+1,lev)))
          CVT % HorizEV4(1:nlongs/2+1,lev) = alpha * CVT % HorizEV4(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV5(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV5(1:nlongs/2+1,lev)))
          CVT % HorizEV5(1:nlongs/2+1,lev) = alpha * CVT % HorizEV5(1:nlongs/2+1,lev)
          alpha                            = SQRT((REAL(nlongs)/2.0) / SUM(CVT % HorizEV6(1:nlongs/2+1,lev) * &
                                                                           CVT % HorizEV6(1:nlongs/2+1,lev)))
          CVT % HorizEV6(1:nlongs/2+1,lev) = alpha * CVT % HorizEV6(1:nlongs/2+1,lev)
        END DO
      END IF

      ! Output
      CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
      PRINT *, 'Outputting horizontal spectra in CVT file: ', TRIM(CVT_filename)
      CALL Write_Covs (CVT_filename, CVT,                     &
                       dims % longs_v(1:nlongs),              &
                       dims % half_levs(1:nlevs),             &
                       output_horiz = .TRUE.)




      PRINT *, 'Computing the vertical transform (reversed MetO transform order)'

      DO loop = 1, Nloops
        PRINT *, '  Finding the mean, loop', loop
        item = 0
        ! Initialise the running total
        CALL Initialise_CVs (ControlVar_mean, .FALSE.)
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Perform the inverse horizontal transform
            CALL U_h_inv(ControlVar_1, ControlVar, CVT, .FALSE.)

            ! Increment the running total
            CALL Add_CVs (ControlVar_mean, ControlVar_1)
          END DO
        END DO

        ! Find the mean for this set of perturbations
        CALL Div_CV_cons (ControlVar_mean, Neffmems_r)

        ! Repeat the reading-in to find the standard deviations
        PRINT *, '  Finding the vertical covs, loop', loop
        item = 0
        DO mem = 1, Nstates
          DO lat = 1, Nlats
            item = item + 1

            ! Read-in this perturbation
            WRITE (paramfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirConParams), '/PertParam_', loop, '_Item', item, '.nc'
            CALL Read_CV (paramfile, ControlVar, 1, CVT, .FALSE., .TRUE.)

            ! Divide by the standard deviations
            CALL U_stddev_inv (ControlVar, CVT)

            ! Perform the inverse horizontal transform
            CALL U_h_inv(ControlVar_1, ControlVar, CVT, .FALSE.)

            ! Remove the mean found above
            CALL Subtract_CVs (ControlVar_1, ControlVar_mean)

            ! Compute contribution to the vertical covs
            CALL CVT_Calibration_vertcovs (CVT, ControlVar_1, zero)  !0.0 means contribute to vert covs

          END DO
        END DO

      END DO

      ! Normalize for vertical covariances
      PRINT *, '  Finding the vert covs'
      CALL CVT_Calibration_vertcovs (CVT, ControlVar, NTotalStates_r)  !ControlVar is not used in this call

      IF (ForceCor) THEN
        ! Due to approximations made, the vertical matrices will not precisely be correlation matrices
        ! Scale to make sure that they are correlation matrices
        DO k = 1, nlongs/2+1
          PRINT *, 'Wavenumber ', k
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode1(1:nlevs,1:nlevs,k))
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode2(1:nlevs,1:nlevs,k))
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode3(1:nlevs,1:nlevs,k))
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode4(1:nlevs,1:nlevs,k))
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode5(1:nlevs,1:nlevs,k))
          CALL Ensure_correlation_matrix (nlevs, CVT % VertMode6(1:nlevs,1:nlevs,k))
        END DO
      END IF

      ! Now compute the vertical eigenvectors and eigenvalues
      PRINT *, '  Finding the eigens'
      CALL VertEigens (CVT)

      ! Output
      CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
      PRINT *, 'Outputting vertical covariances in CVT file: ', TRIM(CVT_filename)
      CALL Write_Covs (CVT_filename, CVT,                     &
                       dims % longs_v(1:nlongs),              &
                       dims % half_levs(1:nlevs),             &
                       output_vert = .TRUE.)


    END SELECT


  END SELECT




END PROGRAM Master_Calibration
