SUBROUTINE ModelObservations ( steps, NLStates, dims, timestep, t0, &
                               Obs, simulate_obs )

!*********************************************************************************
!*                                                                               *
!*  Compute model's version of observations based on pre-run model states        *
!*                                                                               *
!*  steps                  - number of time states (excluding 0)                 *
!*  NLStates               - time sequence of NL model reference states          *
!*  dims                   - dimension data                                      *
!*  timestep               - time step size between states (seconds)             *
!*  t0                     - time corresponding to time step 0 (seconds)         *
!*  Obs                    - pointer to the start of the observation linked list *
!*  simulate_obs           - true if want to compute simulated obs (noise added) *
!*                                                                               *
!*   R. Bannister, 1.4da 02-02-2018                                              *
!*                                                                               *
!*********************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  Obs_type


IMPLICIT NONE

! Parameters
!-----------
INTEGER,                 INTENT(IN) :: steps
TYPE(ABC_type),          INTENT(IN) :: NLStates(0:steps)
TYPE(dims_type),         INTENT(IN) :: dims
REAL(ZREAL8),            INTENT(IN) :: timestep
INTEGER,                 INTENT(IN) :: t0
TYPE(Obs_type), POINTER, INTENT(IN) :: Obs
LOGICAL,                 INTENT(IN) :: simulate_obs

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob
INTEGER                             :: times(0:steps)
INTEGER                             :: t
INTEGER                             :: tsteps_needed4_thisob
REAL(ZREAL8)                        :: prelim_yval, prelim_yval_u, prelim_yval_v, prelim_yval_w

! Functions
!----------
INTEGER                             :: FindLowerIndex
REAL(ZREAL8)                        :: Interpolate3D
REAL(ZREAL8)                        :: GAUSS


! Set-up the times array
DO t=0, steps
  times(t) = t * INT(timestep) + t0
END DO

!PRINT *, 'Inside ModelObservations'
!PRINT *, '------------------------'
!PRINT *, 'time steps ', times(:)

! Loop through observations
thisob => Obs
DO
  IF (ASSOCIATED(thisob)) THEN
    !PRINT *
    !PRINT '(I6,I8,F12.3,F12.3,I3,F12.3)', thisob % batch, thisob % t, thisob % longitude_deg, &
    !               thisob % level_ht, thisob % ob_of_what, thisob % stddev

    ! Find the lower time index
    ! -------------------------
    IF (thisob % tstep_lower == -1) &  ! If not already known
      thisob % tstep_lower = INT(REAL((thisob % t - t0)) / timestep)
    !PRINT *, 'time_step_lower', thisob % tstep_lower

    IF ((thisob % tstep_lower .GE. 0) .AND. (thisob % tstep_lower .LE. steps)) THEN

      ! Determine the number of timesteps needed for this observation
      ! If at beginning or somewhere in the middle of time sequence, then need 2 time steps for linear interp
      ! If at the end, then need only 1
      IF (thisob % tstep_lower .EQ. steps) THEN
        tsteps_needed4_thisob = 1
      ELSE
        tsteps_needed4_thisob = 2
      END IF

      !PRINT *, 'Number of time steps needed for this ob', tsteps_needed4_thisob

      ! Find the lower longitude index
      ! ------------------------------
      IF (thisob % xbox_lower == -1) THEN  ! If not already known
        SELECT CASE (thisob % ob_of_what)
        CASE (1) ! u longitudes
          thisob % xbox_lower    = FindLowerIndex(thisob % longitude_deg, nlongs+2, dims % longs_u(0:nlongs+1))
          thisob % xbox_lower_ws = 0
        CASE (2,3,4,5,6) ! v longitudes
          thisob % xbox_lower    = FindLowerIndex(thisob % longitude_deg, nlongs+2, dims % longs_v(0:nlongs+1))
          thisob % xbox_lower_ws = 0
        CASE (7,8) !Wind speed, so u and v information required
          thisob % xbox_lower    = FindLowerIndex(thisob % longitude_deg, nlongs+2, dims % longs_u(0:nlongs+1))
          thisob % xbox_lower_ws = FindLowerIndex(thisob % longitude_deg, nlongs+2, dims % longs_v(0:nlongs+1))
        END SELECT
      END IF


      IF ((thisob % xbox_lower .GE. 0) .AND. (thisob % xbox_lower_ws .GE. 0)) THEN
        ! Find the lower level index
        ! --------------------------
        IF (thisob % zbox_lower == -1) THEN  ! If not already known
          SELECT CASE (thisob % ob_of_what)
          CASE (1,2,4,7) ! Half levels
            thisob % zbox_lower    = FindLowerIndex(thisob % level_ht, nlevs+2, dims % half_levs(0:nlevs+1))
            thisob % zbox_lower_ws = 0
          CASE (3,5,6) ! Full levels
            thisob % zbox_lower    = FindLowerIndex(thisob % level_ht, nlevs+2, dims % full_levs(0:nlevs+1))
            thisob % zbox_lower_ws = 0
          CASE (8) !Wind speed, so u, v, and w information required
            thisob % zbox_lower    = FindLowerIndex(thisob % level_ht, nlevs+2, dims % half_levs(0:nlevs+1))
            thisob % zbox_lower_ws = FindLowerIndex(thisob % level_ht, nlevs+2, dims % full_levs(0:nlevs+1))
          END SELECT
        END IF

        IF ((thisob % zbox_lower .GE. 0) .AND. (thisob % zbox_lower_ws .GE. 0)) THEN
          ! This observation's model value can be found
          thisob % ob_ok = .TRUE.


          ! Find the observation value (allow non-linear operators based on NLStates data)
          ! ------------------------------------------------------------------------------
          SELECT CASE (thisob % ob_of_what)
          CASE(1)  ! Zonal wind
          !--------------------
            prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
                                         dims % longs_u(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower,                                                       &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         thisob % ob_of_what)
          CASE(2,4)  ! Meridional wind or density
          !--------------------------------------
            prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
                                         dims % longs_v(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower,                                                       &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         thisob % ob_of_what)
          CASE(3,5,6)  ! Vertical wind, buoyancy, or tracer
          !------------------------------------------------
            prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
                                         dims % longs_v(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                         dims % full_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower,                                                       &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         thisob % ob_of_what)
          CASE(7)  ! Horizontal wind speed
          !-------------------------------
            ! Need two pieces of information
            prelim_yval_u = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                           NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                           dims % longs_u(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                           dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                           times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                           thisob % xbox_lower,                                                       &
                                           thisob % zbox_lower,                                                       &
                                           thisob % longitude_deg,                                                    &
                                           thisob % level_ht,                                                         &
                                           thisob % t,                                                                &
                                           1)
            prelim_yval_v = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                           NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                           dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                           dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                           times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                           thisob % xbox_lower_ws,                                                    &
                                           thisob % zbox_lower,                                                       &
                                           thisob % longitude_deg,                                                    &
                                           thisob % level_ht,                                                         &
                                           thisob % t,                                                                &
                                           2)
          CASE(8)  ! Total wind speed
          !--------------------------
            ! Need three pieces of information
            prelim_yval_u = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                           NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                           dims % longs_u(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                           dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                           times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                           thisob % xbox_lower,                                                       &
                                           thisob % zbox_lower,                                                       &
                                           thisob % longitude_deg,                                                    &
                                           thisob % level_ht,                                                         &
                                           thisob % t,                                                                &
                                           1)
            prelim_yval_v = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                           NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                           dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                           dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                           times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                           thisob % xbox_lower_ws,                                                    &
                                           thisob % zbox_lower,                                                       &
                                           thisob % longitude_deg,                                                    &
                                           thisob % level_ht,                                                         &
                                           thisob % t,                                                                &
                                           2)
            prelim_yval_w = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                           NLStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                           dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                           dims % full_levs(thisob % zbox_lower_ws:thisob % zbox_lower_ws+1),         &
                                           times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                           thisob % xbox_lower_ws,                                                    &
                                           thisob % zbox_lower_ws,                                                    &
                                           thisob % longitude_deg,                                                    &
                                           thisob % level_ht,                                                         &
                                           thisob % t,                                                                &
                                           3)
            prelim_yval   = SQRT(prelim_yval_u*prelim_yval_u +    &
                                 prelim_yval_v*prelim_yval_v +    &
                                 prelim_yval_w*prelim_yval_w)
          END SELECT


          !PRINT *, 'prelim_yval', prelim_yval

          ! For what purpose have we generated this model observation?
          IF (simulate_obs) THEN
            ! A simulated observation - add noise to the observation
            thisob % y_true_known = .TRUE.
            thisob % y_true       = prelim_yval
            thisob % y            = prelim_yval + GAUSS(thisob % stddev)
            !PRINT *, 'observation', thisob % y
            !PRINT *
          ELSE
            ! A model observation
            thisob % y_ref        = prelim_yval
            ! Compute the difference between ob and ref
            thisob % d            = thisob % y - thisob % y_ref
            !PRINT *, thisob % y, thisob % y_ref, thisob % d
          END IF


        END IF
      END IF
    END IF

    thisob => thisob % next
  ELSE
    EXIT
  END IF
END DO

END SUBROUTINE ModelObservations
