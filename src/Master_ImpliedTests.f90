PROGRAM Master_ImpliedTests

!*****************************************************
!*   Code to compute test implied cov of vert or     *
!*   horiz transforms only.                          *
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
TYPE(ABC_type)           :: LS
TYPE(CV_type)            :: CV1, CV2, CVinterim, CVinterim1, CVinterim2
TYPE(CVT_type)           :: CVT
INTEGER                  :: point

CHARACTER(LEN=320)       :: LS_filename, CVT_filename, ImplCov_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_ImpliedTests'
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


  ! ============================================================================================
  ! HORIZONTAL TRANSFORM ONLY
  ! ============================================================================================

  IF (ImplCov_npoints > 0) THEN
    PRINT *, '===== Running some point-by-point covariances for the horizontal transform only ====='

    DO point = 1, ImplCov_npoints

      CALL Initialise_CVs (CV1, .FALSE.)
      CALL Initialise_CVs (CV2, .FALSE.)
      CALL Initialise_CVs (CVinterim, .FALSE.)

      ! Put delta function in all fields
      CV1 % v1(longindex(point), levindex(point)) = 1.0
      CV1 % v2(longindex(point), levindex(point)) = 1.0
      CV1 % v3(longindex(point), levindex(point)) = 1.0
      CV1 % v4(longindex(point), levindex(point)) = 1.0
      CV1 % v5(longindex(point), levindex(point)) = 1.0
      CV1 % v6(longindex(point), levindex(point)) = 1.0

      ! Operate with Uh transpose
      CALL U_h_adj (CVinterim, CV1, CVT)
      ! Operate with Uh
      CALL U_h (CVinterim, CV2, CVT)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Horiz_Point_', point, '.nc'
      CALL Write_CV (ImplCov_filename,         &
                     CV2,                      &
                     1,                        &  ! Assume long/lev grid
                     CVT,                      &
                     dims % longs_v(1:nlongs), &
                     dims % full_levs(1:nlevs))

    END DO
  END IF


  ! ============================================================================================
  ! VERTICAL TRANSFORM ONLY
  ! ============================================================================================

  IF (ImplCov_npoints > 0) THEN
    PRINT *, '===== Running some point-by-point covariances for the vertical transform only ====='

    DO point = 1, ImplCov_npoints

      CALL Initialise_CVs (CV1, .FALSE.)
      CALL Initialise_CVs (CV2, .FALSE.)
      CALL Initialise_CVs (CVinterim, .FALSE.)

      ! Put delta function in all fields
      CV1 % v1(longindex(point), levindex(point)) = 1.0
      CV1 % v2(longindex(point), levindex(point)) = 1.0
      CV1 % v3(longindex(point), levindex(point)) = 1.0
      CV1 % v4(longindex(point), levindex(point)) = 1.0
      CV1 % v5(longindex(point), levindex(point)) = 1.0
      CV1 % v6(longindex(point), levindex(point)) = 1.0

      ! Operate with Uv transpose
      CALL U_v_adj (CVinterim, CV1, CVT)
      ! Operate with Uh
      CALL U_v (CVinterim, CV2, CVT)
      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Vert_Point_', point, '.nc'
      CALL Write_CV (ImplCov_filename,         &
                     CV2,                      &
                     1,                        &  ! Assume long/lev grid
                     CVT,                      &
                     dims % longs_v(1:nlongs), &
                     dims % full_levs(1:nlevs))

    END DO
  END IF


  ! ============================================================================================
  ! HORIZONTAL AND VERTICAL TRANSFORMS ONLY
  ! ============================================================================================

  IF (ImplCov_npoints > 0) THEN
    PRINT *, '===== Running some point-by-point covariances for the horizontal and vertical transforms only ====='

    DO point = 1, ImplCov_npoints

      CALL Initialise_CVs (CV1, .FALSE.)
      CALL Initialise_CVs (CV2, .FALSE.)
      CALL Initialise_CVs (CVinterim1, .FALSE.)
      CALL Initialise_CVs (CVinterim2, .FALSE.)

      ! Put delta function in all fields
      CV1 % v1(longindex(point), levindex(point)) = 1.0
      CV1 % v2(longindex(point), levindex(point)) = 1.0
      CV1 % v3(longindex(point), levindex(point)) = 1.0
      CV1 % v4(longindex(point), levindex(point)) = 1.0
      CV1 % v5(longindex(point), levindex(point)) = 1.0
      CV1 % v6(longindex(point), levindex(point)) = 1.0


      SELECT CASE (CVT % CVT_order)

      CASE(1) ! The original Met Office transform order
              ! ---------------------------------------
        CALL U_v_adj (CVinterim1, CV1, CVT)
        CALL U_h_adj (CVinterim2, CVinterim1, CVT)

        CALL U_h (CVinterim2, CVinterim1, CVT)
        CALL U_v (CVinterim1, CV2, CVT)


      CASE(2) ! The reversed horiz/vert transform order
              ! ---------------------------------------
        CALL U_h_adj (CVinterim1, CV1, CVT)
        CALL U_v_adj (CVinterim2, CVinterim1, CVT)

        CALL U_v (CVinterim2, CVinterim1, CVT)
        CALL U_h (CVinterim1, CV2, CVT)

      END SELECT

      ! Output the result
      WRITE (ImplCov_filename, '(A,A,I0.3,A)') TRIM(datadirImpliedCov), '/Spatial_Point_', point, '.nc'
      CALL Write_CV (ImplCov_filename,         &
                     CV2,                      &
                     1,                        &  ! Assume long/lev grid
                     CVT,                      &
                     dims % longs_v(1:nlongs), &
                     dims % full_levs(1:nlevs))

    END DO
  END IF



END PROGRAM Master_ImpliedTests
