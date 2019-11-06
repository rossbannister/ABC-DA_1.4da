SUBROUTINE DeAllocate_Obs(Obs)

! Deallocate a linked list of observations

USE DefConsTypes, ONLY :     &
    Obs_type


IMPLICIT NONE

TYPE(Obs_type), POINTER, INTENT(INOUT) :: Obs

TYPE(Obs_type), POINTER                :: thisob, prevob

IF (ASSOCIATED(Obs)) THEN
  thisob => Obs % next
  DO
    IF (ASSOCIATED(thisob)) THEN
      prevob => thisob
      thisob => thisob % next
      DEALLOCATE(prevob)
    ELSE
      EXIT
    END IF
  END DO
  DEALLOCATE(Obs)
END IF

END SUBROUTINE DeAllocate_Obs
