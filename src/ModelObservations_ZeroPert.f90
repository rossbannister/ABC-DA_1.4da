SUBROUTINE ModelObservations_ZeroPert (Obs)

!*********************************************************************************
!*                                                                               *
!*  Set appropriate parts of observations structure for zero pert state          *
!*  Used when new reference state is used                                        *
!*                                                                               *
!*  Obs                    - pointer to the start of the observation linked list *
!*                                                                               *
!*   R. Bannister, 1.4da 17-03-2018                                              *
!*                                                                               *
!*********************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  Obs_type


IMPLICIT NONE

! Parameters
!-----------
TYPE(Obs_type), POINTER, INTENT(IN) :: Obs

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob



!PRINT *, 'Inside ModelObservations_ZeroPert'
!PRINT *, '---------------------------------'
!PRINT *, 'time steps ', times(:)


! Loop through observations
thisob => Obs
DO
  IF (ASSOCIATED(thisob)) THEN
    IF (thisob % ob_ok) THEN

      thisob % deltay_m     = 0.0
      thisob % hxmy         = -1.0 * thisob % d
      thisob % deltay_m_hat = thisob % hxmy / thisob % variance

    END IF

    thisob => thisob % next
  ELSE
    EXIT
  END IF
END DO

END SUBROUTINE ModelObservations_ZeroPert
