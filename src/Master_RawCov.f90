PROGRAM Master_RawCov

!*****************************************************
!*   Code to compute a selection of raw              *
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
    datadirRawCov,               &
    datadirABCperts,             &
    ImplCov_npoints,             &
    longindex,                   &
    levindex,                    &
    Npointsmax,                  &
    NEnsmax, NEns,               &
    NEnsMems,                    &
    NNMCmax, NNMC,               &
    Nlats


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: ABC_member, Cov
REAL(ZREAL8)             :: Neffmems_inv, source
INTEGER                  :: x_source, z_source
INTEGER                  :: point, var, ens, mem, lat, item
CHARACTER(LEN=6)         :: varname(1:6)
CHARACTER(LEN=320)       :: ABCfile, RawCov_filename


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_RawCov'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions

  IF ((ImplCov_npoints < 1) .OR. (ImplCov_npoints > Npointsmax)) THEN
    PRINT *, 'ImplCov_npoints must be between 1 and ', Npointsmax
    PRINT *, 'Either run with different ImplCov_npoints, or recompile changing Npointsmax'
    STOP
  END IF

  varname(1) = 'u'
  varname(2) = 'v'
  varname(3) = 'w'
  varname(4) = 'r'
  varname(5) = 'b'
  varname(6) = 'tracer'


  Neffmems_inv = 1.0 / REAL(NEns * NEnsMems * Nlats)

  IF (ImplCov_npoints > 0) THEN
    PRINT *, '===== Running some point-by-point covariances ====='

    ! Loop over source points
    DO point = 1, ImplCov_npoints

      x_source = longindex(point)
      z_source = levindex(point)

      ! Repeat for each source field
      DO var = 1, 6
        CALL Initialise_model_vars (Cov, .FALSE.)
        IF (NEns > 0) THEN
          DO ens = 1, NEns
            item = 0
            DO lat = 1, Nlats
              DO mem = 1, NEnsMems
                item = item + 1
                ! Read-in this ensemble member
                WRITE (ABCfile, '(A,A,I0.3,A,I0.3,A)') TRIM(datadirABCperts), '/PertABC_Ens', ens, '_Item', item, '.nc'
                PRINT *, 'Reading file ', TRIM(ABCfile)
                CALL Read_state_2d (ABCfile, ABC_member, dims, 1)

                IF ((point ==1).AND.(var == 1).AND.(ens == 1).AND.(mem == 1)) THEN
                  ! Set-up model grid with info in file read-in
                  CALL Set_ht_dep_cons (dims)
                END IF

                ! Source point
                SELECT CASE (var)
                CASE (1)
                  source = ABC_member % u(x_source,z_source)
                CASE (2)
                  source = ABC_member % v(x_source,z_source)
                CASE (3)
                  source = ABC_member % w(x_source,z_source)
                CASE (4)
                  source = ABC_member % r(x_source,z_source)
                CASE (5)
                  source = ABC_member % b(x_source,z_source)
                CASE (6)
                  source = ABC_member % tracer(x_source,z_source)
                END SELECT
                Cov % u(1:nlongs,1:nlevs)      = Cov % u(1:nlongs,1:nlevs)      + source * ABC_member % u(1:nlongs,1:nlevs)
                Cov % v(1:nlongs,1:nlevs)      = Cov % v(1:nlongs,1:nlevs)      + source * ABC_member % v(1:nlongs,1:nlevs)
                Cov % w(1:nlongs,1:nlevs)      = Cov % w(1:nlongs,1:nlevs)      + source * ABC_member % w(1:nlongs,1:nlevs)
                Cov % r(1:nlongs,1:nlevs)      = Cov % r(1:nlongs,1:nlevs)      + source * ABC_member % r(1:nlongs,1:nlevs)
                Cov % b(1:nlongs,1:nlevs)      = Cov % b(1:nlongs,1:nlevs)      + source * ABC_member % b(1:nlongs,1:nlevs)
                Cov % tracer(1:nlongs,1:nlevs) = Cov % tracer(1:nlongs,1:nlevs) + source * ABC_member % tracer(1:nlongs,1:nlevs)
              END DO
            END DO
          END DO

          ! Normalize
          Cov % u(1:nlongs,1:nlevs)      = Cov % u(1:nlongs,1:nlevs)      * Neffmems_inv
          Cov % v(1:nlongs,1:nlevs)      = Cov % v(1:nlongs,1:nlevs)      * Neffmems_inv
          Cov % w(1:nlongs,1:nlevs)      = Cov % w(1:nlongs,1:nlevs)      * Neffmems_inv
          Cov % r(1:nlongs,1:nlevs)      = Cov % r(1:nlongs,1:nlevs)      * Neffmems_inv
          Cov % b(1:nlongs,1:nlevs)      = Cov % b(1:nlongs,1:nlevs)      * Neffmems_inv
          Cov % tracer(1:nlongs,1:nlevs) = Cov % tracer(1:nlongs,1:nlevs) * Neffmems_inv

          ! Output these covariance fields
          WRITE (RawCov_filename, '(A,A,I0.3,A,A,A)') TRIM(datadirRawCov), '/Point_', point, &
                 '_delta', TRIM(varname(var)), '.nc'
          CALL Write_state_2d (RawCov_filename, Cov, dims, 1, 0, 0)


        END IF
      END DO

    END DO
  END IF

END PROGRAM Master_RawCov
