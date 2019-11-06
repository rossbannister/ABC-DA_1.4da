SUBROUTINE ModelObservations_linear ( steps, NLStates, PertStates, dims, timestep, t0, &
                                      Obs )

!*********************************************************************************
!*                                                                               *
!*  Compute model's version of observations based on pre-run model states        *
!*  (linear perturbations)                                                       *
!*                                                                               *
!*  steps                  - number of time states (excluding 0)                 *
!*  NLStates               - time sequence of NL model reference states          *
!*  PertStates             - time sequence of pert states                        *
!*  dims                   - dimension data                                      *
!*  timestep               - time step size between states (seconds)             *
!*  t0                     - time corresponding to time step 0 (seconds)         *
!*  Obs                    - pointer to the start of the observation linked list *
!*                                                                               *
!*   R. Bannister, 1.4da 17-03-2018                                              *
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
TYPE(ABC_type),          INTENT(IN) :: PertStates(0:steps)
TYPE(dims_type),         INTENT(IN) :: dims
REAL(ZREAL8),            INTENT(IN) :: timestep
INTEGER,                 INTENT(IN) :: t0
TYPE(Obs_type), POINTER, INTENT(IN) :: Obs

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob
INTEGER                             :: times(0:steps)
INTEGER                             :: t
INTEGER                             :: tsteps_needed4_thisob
REAL(ZREAL8)                        :: prelim_yval, prelim_yval_u, prelim_yval_v, prelim_yval_w
REAL(ZREAL8)                        :: d_prelim_yval, d_prelim_yval_u, d_prelim_yval_v, d_prelim_yval_w

! Functions
!----------
INTEGER                             :: FindLowerIndex
REAL(ZREAL8)                        :: Interpolate3D
REAL(ZREAL8)                        :: GAUSS


! Set-up the times array
DO t=0, steps
  times(t) = t * INT(timestep) + t0
END DO

!PRINT *, 'Inside ModelObservations_linear'
!PRINT *, '-------------------------------'
!PRINT *, 'time steps ', times(:)


! Loop through observations
thisob => Obs
DO
  IF (ASSOCIATED(thisob)) THEN
    IF (thisob % ob_ok) THEN
      ! This indicates that the tstep_lower, xbox_lower, and zbox_lower are already computed

      ! Determine the number of timesteps needed for this observation
      ! If at beginning or somewhere in the middle of time sequence, then need 2 time steps for linear interp
      ! If at the end, then need only 1
      IF (thisob % tstep_lower .EQ. steps) THEN
        tsteps_needed4_thisob = 1
      ELSE
        tsteps_needed4_thisob = 2
      END IF

      ! Find the observation perturbations based on model perturbations
      ! ---------------------------------------------------------------
      SELECT CASE (thisob % ob_of_what)
      CASE(1)  ! Zonal wind
      !--------------------
        d_prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                       PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
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
        d_prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                       PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
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
        d_prelim_yval = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                       PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1), &
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
        d_prelim_yval_u = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                         dims % longs_u(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower,                                                       &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         1)
        d_prelim_yval_v = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                         dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower_ws,                                                    &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         2)
        d_prelim_yval   = ( prelim_yval_u * d_prelim_yval_u +    &
                            prelim_yval_v * d_prelim_yval_v ) / thisob % y_ref

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
        d_prelim_yval_u = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                         dims % longs_u(thisob % xbox_lower:thisob % xbox_lower+1),                 &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower,                                                       &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         1)
        d_prelim_yval_v = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                         dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                         dims % half_levs(thisob % zbox_lower:thisob % zbox_lower+1),               &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower_ws,                                                    &
                                         thisob % zbox_lower,                                                       &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         2)
        d_prelim_yval_w = Interpolate3D (tsteps_needed4_thisob,                                                     &
                                         PertStates(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),&
                                         dims % longs_v(thisob % xbox_lower_ws:thisob % xbox_lower_ws+1),           &
                                         dims % full_levs(thisob % zbox_lower_ws:thisob % zbox_lower_ws+1),         &
                                         times(thisob % tstep_lower:thisob % tstep_lower+tsteps_needed4_thisob-1),  &
                                         thisob % xbox_lower_ws,                                                    &
                                         thisob % zbox_lower_ws,                                                    &
                                         thisob % longitude_deg,                                                    &
                                         thisob % level_ht,                                                         &
                                         thisob % t,                                                                &
                                         3)
        d_prelim_yval   = ( prelim_yval_u * d_prelim_yval_u +    &
                            prelim_yval_v * d_prelim_yval_v +    &
                            prelim_yval_w * d_prelim_yval_w ) / thisob % y_ref
      END SELECT


      thisob % deltay_m     = d_prelim_yval
      thisob % hxmy         = thisob % deltay_m - thisob % d
      thisob % deltay_m_hat = thisob % hxmy / thisob % variance

    END IF

    thisob => thisob % next
  ELSE
    EXIT
  END IF
END DO

END SUBROUTINE ModelObservations_linear
