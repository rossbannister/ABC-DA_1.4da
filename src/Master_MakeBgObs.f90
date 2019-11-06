PROGRAM Master_MakeBgObs

!*****************************************************
!*   Code to make synthetic background and obs       *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    dims_type,                   &
    ABC_type,                    &
    CV_type,                     &
    CVT_type,                    &
    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    init_ABC_file,               &
    output_ABC_file,             &
    ntimesteps,                  &
    ndumps,                      &
    diagnostics_file,            &
    Obs_type,                    &
    ObsSpec_type,                &
    Generate_mode,               &
    ObsSpec,                     &
    datadir_ObsSpec,             &
    ObsSpec_file,                &
    datadir_Obs,                 &
    datadir_Bg,                  &
    Bg_file,                     &
    Obs_file,                    &
    datadirABC_in,               &
    dt,                          &
    dt_da,                       &
    runlength,                   &
    zero,                        &
    t0,                          &
    Pert_file

IMPLICIT NONE

INCLUDE "DeAllocate_Obs.interface"
INCLUDE "Boundaries_CV.interface"
INCLUDE "ModelObservations.interface"
INCLUDE "Write_Obs.interface"

! Declare variables
!==========================
TYPE(dims_type)             :: dims
TYPE(ABC_type), ALLOCATABLE :: ABC(:)
TYPE(CV_type)               :: CV_data
TYPE(CVT_type)              :: CVT

CHARACTER(LEN=320)          :: GenerateObs_filename, Obs_filename, ABC_init_filename
CHARACTER(LEN=320)          :: Truth_filename, CVT_filename, Pert_filename, Bg_filename
INTEGER                     :: batch, obcount, x, z, t
REAL(ZREAL8)                :: longitude, height
CHARACTER(LEN=320)          :: blank
INTEGER                     :: version, IOstatus, linecount
TYPE(Obs_type), POINTER     :: Observations, thisob
LOGICAL                     :: FirstObRead
INTEGER                     :: maxtime, timestepsreq, DAtimestepsreq, model_steps_per_da_step


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_MakeBgObs'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions

  SELECT CASE (Generate_mode)

  CASE (1)
    ! Generate file that specifies obs positions/times/etc (not the obs themselves)
    ! -----------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------

    ! Open the output file for the obs specification file
    GenerateObs_filename = TRIM(datadir_ObsSpec) // '/' // TRIM(ObsSpec_file)
    PRINT *, 'Generating observations in file ', TRIM(GenerateObs_filename)

    OPEN (13, file=GenerateObs_filename)
    WRITE (13,'(A)') 'Observation specification file for ABC model'
    WRITE (13,'(A18,I3)') 'Format version  : ', 1
    WRITE (13,'(A)') '-----------------------------------------'
    WRITE (13,'(A18,I4)') 'Ref year        : ', ObsSpec % year0
    WRITE (13,'(A18,I2)') 'Ref month       : ', ObsSpec % month0
    WRITE (13,'(A18,I2)') 'Ref day         : ', ObsSpec % day0
    WRITE (13,'(A18,I2)') 'Ref hour        : ', ObsSpec % hour0
    WRITE (13,'(A18,I2)') 'Ref minute      : ', ObsSpec % min0
    WRITE (13,'(A18,I2)') 'Ref second      : ', ObsSpec % sec0

    obcount = 0

    DO batch = 1, ObsSpec % NumBatches
      DO x = 1, ObsSpec % NumObs_long(batch)
        longitude = ObsSpec % long_min(batch) +                                   &
                    (ObsSpec % long_max(batch) - ObsSpec % long_min(batch)) *     &
                    REAL(x-1) / REAL(ObsSpec % NumObs_long(batch)-1)
        DO z = 1, ObsSpec % NumObs_height(batch)
          height = ObsSpec % height_min(batch) +                                  &
                    (ObsSpec % height_max(batch) - ObsSpec % height_min(batch)) * &
                    REAL(z-1) / REAL(ObsSpec % NumObs_height(batch)-1)

          obcount = obcount + 1
          WRITE (13,'(A)') '-----------------------------------------'
          WRITE (13,'(A18,I10)')   'Observation No  : ', obcount
          WRITE (13,'(A18,I6)')    'Batch ID        : ', ObsSpec % batch(batch)
          WRITE (13,'(A18,I8)')    'Time of obs (s) : ', ObsSpec % seconds(batch)
          WRITE (13,'(A18,F12.3)') 'Longitude (deg) : ', longitude
          WRITE (13,'(A18,F12.3)') 'Height (m)      : ', height
          WRITE (13,'(A18,I3)')    'Observation of  : ', ObsSpec % ob_of_what(batch)
          WRITE (13,'(A18,F12.8)') 'Err stddev      : ', ObsSpec % stddev(batch)
        END DO
      END DO
    END DO


    CLOSE (13)


  CASE (2)
    ! Generate synthetic observations according to specified file
    ! -----------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------

    ! File containing the obs specifications
    GenerateObs_filename = TRIM(datadir_ObsSpec) // '/' // TRIM(ObsSpec_file)
    PRINT *, 'Obs Specification file ', TRIM(GenerateObs_filename)

    Obs_filename = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
    PRINT *, 'Observations file ', TRIM(Obs_filename)

    ! Initialise the observation structure (a linked list)
    NULLIFY(Observations)

    OPEN (13, file=GenerateObs_filename)
    READ (13,'(A)') blank
    READ (13,'(A18,I3)') blank, version
    PRINT *, 'File is of version ', version
    READ (13,*) blank
    READ (13,'(A18,I4)') blank, ObsSpec % year0
    READ (13,'(A18,I2)') blank, ObsSpec % month0
    READ (13,'(A18,I2)') blank, ObsSpec % day0
    READ (13,'(A18,I2)') blank, ObsSpec % hour0
    READ (13,'(A18,I2)') blank, ObsSpec % min0
    READ (13,'(A18,I2)') blank, ObsSpec % sec0
    PRINT '(A,I5,I3,I3,I3,I3,I3)', 'Ref time (year, month, day, hour, min, sec) ', &
             ObsSpec % year0, ObsSpec % month0, ObsSpec % day0,                    &
             ObsSpec % hour0, ObsSpec % min0, ObsSpec % sec0
    linecount   = 9
    FirstObRead = .FALSE.
    maxtime     = t0
    DO
      linecount = linecount + 1
      READ(13,'(A)',IOSTAT=IOstatus)
      IF (IOstatus >0) THEN
        PRINT*, 'Error reading observation specification file at line ', linecount
        STOP
      ELSE
        IF (IOstatus < 0) EXIT    ! Exit loop
        ! Make a record for this observation and read-in details
        IF (.NOT.FirstObRead) THEN
          ALLOCATE (Observations)
          thisob      => Observations
          FirstObRead = .TRUE.
        ELSE
          ALLOCATE (thisob % next)
          thisob => thisob % next
        END IF
        READ (13,'(A18,I10)')   blank, thisob % obnumber_thisfile
        READ (13,'(A18,I6)')    blank, thisob % batch
        READ (13,'(A18,I8)')    blank, thisob % t
        READ (13,'(A18,F12.3)') blank, thisob % longitude_deg
        READ (13,'(A18,F12.3)') blank, thisob % level_ht
        READ (13,'(A18,I3)')    blank, thisob % ob_of_what
        READ (13,'(A18,F12.8)') blank, thisob % stddev
        PRINT '(I10,I6,I8,F12.3,F12.3,I3,F12.3)', thisob % obnumber_thisfile, thisob % batch, thisob % t, thisob % longitude_deg, &
                 thisob % level_ht, thisob % ob_of_what, thisob % stddev
        ! Initialize the rest of the structure
        thisob % ob_ok         = .FALSE.
        !thisob % year          = 0
        !thisob % month         = 0
        !thisob % day           = 0
        !thisob % hour          = 0
        !thisob % min           = 0
        !thisob % sec           = 0
        thisob % xbox_lower    = -1
        thisob % xbox_lower_ws = -1
        thisob % zbox_lower    = -1
        thisob % zbox_lower_ws = -1
        thisob % tstep_lower   = -1
        thisob % y_true_known  = .FALSE.
        thisob % y_true        = 0.0
        thisob % y             = 0.0
        thisob % variance      = thisob % stddev * thisob % stddev
        thisob % y_ref         = 0.0
        thisob % d             = 0.0
        thisob % deltay_m      = 0.0
        thisob % hxmy          = 0.0
        thisob % deltay_m_hat  = 0.0
        NULLIFY(thisob % next)
        IF (thisob % t > maxtime) maxtime = thisob % t
      END IF
    END DO

    CLOSE (13)


    ! Run the non-linear model to generate the synthetic observations
    ! ---------------------------------------------------------------

    PRINT *, 'The latest observation is at', maxtime
    maxtime = maxtime - t0
    PRINT *, '  ... since t0              ', maxtime
    PRINT *, 'The model time step is', dt
    PRINT *, 'The da time step is', dt_da
    PRINT *, 'The specified run length is ', INT(runlength)

    IF (REAL(maxtime) > runlength) THEN
      PRINT *, 'Error: One or more obs times specified exceed the window length'
      STOP
    END IF

    ! How many DA timesteps?
    DAtimestepsreq = INT(runlength / dt_da + 0.999999)
    PRINT *, 'Number of DA timesteps ', DAtimestepsreq

    ! How many model timesteps?
    timestepsreq   = DAtimestepsreq * dt_da / dt
    PRINT *, 'Number of model timesteps ', timestepsreq

    !Allocate the array of data of DA time steps
    ALLOCATE (ABC(0:DAtimestepsreq))

    ! The truth at the initial time
    ABC_init_filename   = TRIM(datadirABC_in) // '/' // TRIM(init_ABC_file)
    PRINT *, 'The file containing the truth: ', TRIM(ABC_init_filename)

    ! Set state to zero
    CALL Initialise_model_vars (ABC(0), .FALSE.)
    CALL Initialise_dims (dims)

    ! Read in initial conditions
    PRINT*, 'Reading in initial conditions ...'
    CALL Read_state_2d (ABC_init_filename, ABC(0), dims, -1)
    PRINT*, '-- done'

    ! Set some commonly-used constants
    CALL Set_ht_dep_cons (dims)
    PRINT *, 'About to run model for ', timestepsreq, ' time steps'
    CALL ABC_NL_ModelDriver_DA (ABC(0:DAtimestepsreq), dims, timestepsreq, DAtimestepsreq)
    PRINT *, ' -- done'

    ! Output the truth trajectory
    Truth_filename = TRIM(datadir_Obs) // '/' // TRIM(output_ABC_file)
    model_steps_per_da_step = INT(dt_da/dt)
    DO t = 0, DAtimestepsreq
      ! Calculate diagnostics
      CALL Calc_geost (ABC(t))
      CALL Calc_hydro (ABC(t), dims)
      CALL Calc_vert_mom_source (ABC(t), dims)
      CALL Calc_horiz_div(ABC(t), dims)
      CALL Calc_horiz_vort(ABC(t), dims)
      CALL Effective_buoyancy (ABC(t), dims)
      CALL Write_state_2d (Truth_filename, ABC(t), dims, DAtimestepsreq+1, t, model_steps_per_da_step)
    END DO

    ! Generate the model observations
    PRINT *, 'About to generate model observations'
    CALL ModelObservations (DAtimestepsreq, ABC(0:DAtimestepsreq), dims, dt_da, t0, &
                            Observations, .TRUE.)
    PRINT *, ' -- done'


    ! Output the observations
    PRINT *, 'Outputting the observations'
    Obs_filename = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
    PRINT *, 'Observations in file ', TRIM(Obs_file)
    CALL Write_Obs (Obs_filename, INT(runlength), dt, timestepsreq, dt_da, DAtimestepsreq, &
                    Observations)

    ! Deallocate the observation structure
    CALL DeAllocate_Obs (Observations)
    ! Deallocate the model fields
    DEALLOCATE (ABC)


  CASE (3)
    ! Generate a background state for the start of a cycle
    ! -----------------------------------------------------------------------------
    ! -----------------------------------------------------------------------------

    ALLOCATE (ABC(1:3))
    ! 1 is for the truth (read-in)
    ! 2 is for the pert (computed)
    ! 3 is for the generated background (sum of above)

    ! The truth at the initial time
    ABC_init_filename = TRIM(datadirABC_in) // '/' // TRIM(init_ABC_file)
    PRINT *, 'The file containing the truth: ', TRIM(ABC_init_filename)

    ! Set state to zero
    CALL Initialise_model_vars (ABC(1), .FALSE.)
    CALL Initialise_dims (dims)

    ! Read in initial conditions
    PRINT*, 'Reading in initial conditions (truth) ...'
    CALL Read_state_2d (ABC_init_filename, ABC(1), dims, -1)
    PRINT*, '-- done'

    ! Set some commonly-used constants
    CALL Set_ht_dep_cons (dims)

    ! The CVT filename
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    PRINT *, 'The CVT file: ', TRIM(CVT_filename)

    PRINT*, 'Reading in CVT data ...'
    CALL Read_Covs (CVT_filename, CVT,              &
                    .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
    PRINT *, '-- done'

    ! Put some IID random numbers (variance 1) in the control variable elements
    CALL Initialise_CVs (CV_data, .TRUE.)

    ! Act with the forward CVT to generate the perturbation
    CALL U_trans (ABC(1), CV_data, ABC(2), CVT, dims)
    ! Compute diagnostics for this pert in model space
    CALL Calc_hydro (ABC(2), dims)
    CALL Calc_geost (ABC(2))
    CALL Calc_vert_mom_source(ABC(2), dims)
    CALL Calc_horiz_div(ABC(2), dims)
    CALL Calc_horiz_vort(ABC(2), dims)
    CALL Effective_buoyancy (ABC(2), dims)


    ! Add the background increment to the truth to generate a background
    ABC(3) = ABC(1)
    CALL Add_model_vars (ABC(3), ABC(2), .TRUE.)
    ! Compute diagnostics for the generated background state
    CALL Calc_hydro (ABC(3), dims)
    CALL Calc_geost (ABC(3))
    CALL Calc_vert_mom_source(ABC(3), dims)
    CALL Calc_horiz_div(ABC(3), dims)
    CALL Calc_horiz_vort(ABC(3), dims)
    CALL Effective_buoyancy (ABC(3), dims)


    ! Output the perturbation and background states
    ! The pert filename
    Pert_filename = TRIM(datadir_Bg) // '/' // TRIM(Pert_file)
    PRINT *, 'The bg pert output file: ', TRIM(Pert_filename)
    CALL Write_state_2d (Pert_filename, ABC(2), dims, 1, 0, 0)

    ! The background filename
    Bg_filename = TRIM(datadir_Bg) // '/' // TRIM(Bg_file)
    PRINT *, 'The bg pert output file: ', TRIM(Bg_filename)
    CALL Write_state_2d (Bg_filename, ABC(3), dims, 1, 0, 0)

    DEALLOCATE (ABC)

  END SELECT


END PROGRAM Master_MakeBgObs
