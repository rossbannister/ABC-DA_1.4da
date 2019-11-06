SUBROUTINE Boundaries_CV (state, set_1, set_2, set_3, set_4, set_5, set_6)

!************************************
!* Subroutine to apply boundary     *
!* conditions to control params     *
!*                                  *
!* Input flags to choose which      *
!* variables to apply boundary      *
!* conditions to                    *
!*                                  *
!* Set no flags to do them all      *
!************************************

USE DefConsTypes, ONLY :     &
    CV_type,                 &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CV_type),        INTENT(INOUT) :: state
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_1
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_2
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_3
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_4
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_5
LOGICAL, OPTIONAL,    INTENT(IN)    :: set_6

!Local variables
!---------------
LOGICAL                             :: flag_1, flag_2, flag_3, flag_4, flag_5, flag_6
LOGICAL                             :: some_flags



IF (PRESENT(set_1)) THEN
  flag_1 = set_1
ELSE
  flag_1 = .FALSE.
END IF
IF (PRESENT(set_2)) THEN
  flag_2 = set_2
ELSE
  flag_2 = .FALSE.
END IF
IF (PRESENT(set_3)) THEN
  flag_3 = set_3
ELSE
  flag_3 = .FALSE.
END IF
IF (PRESENT(set_4)) THEN
  flag_4 = set_4
ELSE
  flag_4 = .FALSE.
END IF
IF (PRESENT(set_5)) THEN
  flag_5 = set_5
ELSE
  flag_5 = .FALSE.
END IF
IF (PRESENT(set_6)) THEN
  flag_6 = set_6
ELSE
  flag_6 = .FALSE.
END IF


some_flags = PRESENT(set_1) .OR. PRESENT(set_2) .OR. PRESENT(set_3) .OR. &
             PRESENT(set_4) .OR. PRESENT(set_5) .OR. PRESENT(set_6)

IF (.NOT.some_flags) THEN
  ! No flags set.  This is shorthand for all flags
  flag_1      = .TRUE.
  flag_2      = .TRUE.
  flag_3      = .TRUE.
  flag_4      = .TRUE.
  flag_5      = .TRUE.
  flag_6      = .TRUE.
END IF


IF (flag_1) THEN   !psi
  ! Horizontal boundaries
  state % v1(0, 1:nlevs)          = state % v1(nlongs, 1:nlevs)
  state % v1(nlongs+1, 1:nlevs)   = state % v1(1, 1:nlevs)
  ! Vertical boundaries
  state % v1(0:nlongs+1,0)        = - 1.0 * state % v1(0:nlongs+1,1)
  state % v1(0:nlongs+1,nlevs+1)  = state % v1(0:nlongs+1, nlevs)
ENDIF

IF (flag_2) THEN   !chi
  ! Horizontal boundaries
  state % v2(0, 1:nlevs)          = state % v2(nlongs, 1:nlevs)
  state % v2(nlongs+1, 1:nlevs)   = state % v2(1, 1:nlevs)
  ! Vertical boundaries
  state % v2(0:nlongs+1,0)        = - 1.0 * state % v2(0:nlongs+1,1)
  state % v2(0:nlongs+1,nlevs+1)  = state % v2(0:nlongs+1, nlevs)
ENDIF

IF (flag_3) THEN   !r-like
  ! Horizontal boundaries
  state % v3(0, 1:nlevs)          = state % v3(nlongs, 1:nlevs)
  state % v3(nlongs+1, 1:nlevs)   = state % v3(1, 1:nlevs)
  ! Vertical boundaries
  state % v3(0:nlongs+1,0)        = state % v3(0:nlongs+1,1)
  state % v3(0:nlongs+1,nlevs+1)  = state % v3(0:nlongs+1,nlevs)
ENDIF

IF (flag_4) THEN   !b-like
  ! Horizontal boundaries
  state % v4(0, 1:nlevs)          = state % v4(nlongs, 1:nlevs)
  state % v4(nlongs+1, 1:nlevs)   = state % v4(1, 1:nlevs)
  ! Vertical boundaries
  state % v4(0:nlongs+1,0)        = 0.0
  state % v4(0:nlongs+1,nlevs)    = 0.0
  state % v4(0:nlongs+1,nlevs+1)  = 0.0
ENDIF

IF (flag_5) THEN   !w-like
  ! Horizontal boundaries
  state % v5(0, 1:nlevs)          = state % v5(nlongs, 1:nlevs)
  state % v5(nlongs+1, 1:nlevs)   = state % v5(1, 1:nlevs)
  ! Vertical boundaries
  state % v5(0:nlongs+1,0)        = 0.0
  state % v5(0:nlongs+1,nlevs)    = 0.0
  state % v5(0:nlongs+1,nlevs+1)  = 0.0
ENDIF

IF (flag_6) THEN   !tracer-like
  ! Horizontal boundaries
  state % v6(0, 1:nlevs)          = state % v6(nlongs, 1:nlevs)
  state % v6(nlongs+1, 1:nlevs)   = state % v6(1, 1:nlevs)
  ! Vertical boundaries
  state % v6(0:nlongs+1, 0)       = state % v6(0:nlongs+1, 1)
  state % v6(0:nlongs+1, nlevs+1) = state % v6(0:nlongs+1, nlevs)
ENDIF

END SUBROUTINE Boundaries_CV
