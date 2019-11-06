SUBROUTINE Read_Obs (Observations, Obs_filename, &
                     dt, timestepsreq, dt_da, DAtimestepsreq, maxtime)

!*********************************************************************************
!*                                                                               *
!*  Read-in observations                                                         *
!*  (linear perturbations)                                                       *
!*                                                                               *
!*                                                                               *
!*   R. Bannister, 1.4da 28-03-2018                                              *
!*                                                                               *
!*********************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  Obs_type


IMPLICIT NONE

! Parameters
!-----------
TYPE(Obs_type), POINTER, INTENT(INOUT) :: Observations
CHARACTER(LEN=320),      INTENT(IN)    :: Obs_filename
REAL(ZREAL8),            INTENT(OUT)   :: dt
INTEGER,                 INTENT(OUT)   :: timestepsreq
REAL(ZREAL8),            INTENT(OUT)   :: dt_da
INTEGER,                 INTENT(OUT)   :: DAtimestepsreq
INTEGER,                 INTENT(OUT)   :: maxtime


! Local Variables
!----------------
TYPE(Obs_type), POINTER                :: thisob
CHARACTER(LEN=320)                     :: blank
INTEGER                                :: obcount, version, IOstatus
LOGICAL                                :: FirstObRead


  OPEN (13, file=Obs_filename)
  READ (13,'(A)')         blank
  READ (13,'(A18,I3)')    blank, version
  READ (13,'(A18,I8)')    blank, maxtime
  READ (13,'(A18,E12.3)') blank, dt
  READ (13,'(A18,I8)')    blank, timestepsreq
  READ (13,'(A18,E12.3)') blank, dt_da
  READ (13,'(A18,I8)')    blank, DAtimestepsreq
  READ (13,'(A)')         blank

  obcount     = 0
  FirstObRead = .FALSE.

  DO

    READ(13,'(A)',IOSTAT=IOstatus) blank
    IF (IOstatus >0) THEN
      PRINT*, 'Error reading observations file at obs ', obcount
      STOP
    ELSE
      IF (IOstatus < 0) EXIT    ! Exit loop
    END IF


    IF (.NOT.FirstObRead) THEN
      ALLOCATE (Observations)
      thisob      => Observations
      FirstObRead = .TRUE.
    ELSE
      ALLOCATE (thisob % next)
      thisob => thisob % next
    END IF

    obcount = obcount + 1
    READ (13,'(A18,I10)')   blank, thisob % obnumber_thisfile
    READ (13,'(A18,I6)')    blank, thisob % batch
    READ (13,'(A18,I8)')    blank, thisob % t
    READ (13,'(A18,F12.3)') blank, thisob % longitude_deg
    READ (13,'(A18,F12.3)') blank, thisob % level_ht
    READ (13,'(A18,I10)')   blank, thisob % xbox_lower
    READ (13,'(A18,I10)')   blank, thisob % xbox_lower_ws
    READ (13,'(A18,I10)')   blank, thisob % zbox_lower
    READ (13,'(A18,I10)')   blank, thisob % zbox_lower_ws
    READ (13,'(A18,I10)')   blank, thisob % tstep_lower
    READ (13,'(A18,I3)')    blank, thisob % ob_of_what
    READ (13,'(A)')         blank
    READ (13,'(A18,L)')     blank, thisob % y_true_known
    READ (13,'(A18,E14.5)') blank, thisob % y_true
    READ (13,'(A18,E14.5)') blank, thisob % y
    READ (13,'(A18,E14.5)') blank, thisob % stddev
    READ (13,'(A18,E14.5)') blank, thisob % y_ref
    READ (13,'(A18,E14.5)') blank, thisob % d
    READ (13,'(A18,E14.5)') blank, thisob % deltay_m
    READ (13,'(A18,E14.5)') blank, thisob % hxmy
    READ (13,'(A18,E14.5)') blank, thisob % deltay_m_hat
    thisob % ob_ok    = .TRUE.
    thisob % variance = thisob % stddev * thisob % stddev
    NULLIFY(thisob % next)

  END DO

  CLOSE (13)

END SUBROUTINE Read_Obs
