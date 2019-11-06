PROGRAM Master_ImpliedCov

!*****************************************************
!*   Code to compute a selection of implied          *
!*   background error covariances.                   *
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
    datadirImpliedCov,           &
    datadirABCfcs,               &
    LS_file,                     &
    datadirCVT,                  &
    CVT_file,                    &
    ImplCov_npoints,             &
    longindex,                   &
    levindex,                    &
    Npointsmax


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: Model_data1, Model_data2, LS
TYPE(CV_type)            :: CV_data
TYPE(CVT_type)           :: CVT
INTEGER                  :: point

CHARACTER(LEN=320)       :: LS_filename, CVT_filename, ImplCov_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_ImpliedCov'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions

  IF ((ImplCov_npoints < 1) .OR. (ImplCov_npoints > Npointsmax)) THEN
    PRINT *, 'ImplCov_npoints must be between 1 and ', Npointsmax
    PRINT *, 'Either run with different ImplCov_npoints, or recompile changing Npointsmax'
    STOP
  END IF


  LS_filename          = TRIM(datadirABCfcs)   // '/' // TRIM(LS_file)
  CVT_filename         = TRIM(datadirCVT)      // '/' // TRIM(CVT_file)

  ! Read in LS
  PRINT*, 'Reading in linearisation state ...'
  CALL Read_state_2d (LS_filename, LS, dims, 1)
  PRINT*, '-- done'

  ! Set some commonly-used constants
  CALL Set_ht_dep_cons (dims)

  ! Read in CVT data prepared beforehand
  PRINT*, 'Reading in CVT data ...'
  CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
  CALL Read_Covs (CVT_filename, CVT,              &
                  .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
  PRINT *, '-- done'




  IF (ImplCov_npoints > 0) THEN
    PRINT *, '===== Running some point-by-point covariances ====='

    DO point = 1, ImplCov_npoints

      CALL Initialise_model_vars (Model_data1, .FALSE.)

      ! Put delta function in the u field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta u'
      Model_data1 % u(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL Initialise_model_vars (Model_data2, .FALSE.)
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltau.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)


      ! Put delta function in the v field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta v'
      Model_data1 % u(longindex(point), levindex(point)) = 0.0
      Model_data1 % v(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL Initialise_CVs (CV_data, .FALSE.)
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltav.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)


      ! Put delta function in the w field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta w'
      Model_data1 % v(longindex(point), levindex(point)) = 0.0
      Model_data1 % w(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL Initialise_CVs (CV_data, .FALSE.)
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltaw.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)


      ! Put delta function in the r field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta r'
      Model_data1 % w(longindex(point), levindex(point)) = 0.0
      Model_data1 % r(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL Initialise_CVs (CV_data, .FALSE.)
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltar.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)


      ! Put delta function in the b field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta b'
      Model_data1 % r(longindex(point), levindex(point)) = 0.0
      Model_data1 % b(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL Initialise_CVs (CV_data, .FALSE.)
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltab.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)


      ! Put delta function in the tracer field
      ! ---------------------------------
      PRINT *, '  Point', point, ' delta b'
      Model_data1 % b(longindex(point), levindex(point)) = 0.0
      Model_data1 % tracer(longindex(point), levindex(point)) = 1.0
      ! Operate with U transpose
      CALL U_trans_adj (LS, CV_data, Model_data1, CVT, dims)
      ! Operate with U
      CALL Initialise_model_vars (Model_data2, .FALSE.)
      CALL U_trans (LS, CV_data, Model_data2, CVT, dims)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Point_', point, '_deltatracer.nc'
      CALL Write_state_2d (ImplCov_filename, Model_data2, dims, 1, 0, 0)

    END DO
  END IF

END PROGRAM Master_ImpliedCov
