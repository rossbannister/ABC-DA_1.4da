SUBROUTINE Write_Obs ( Filename, maxtime, dt, timestepsreq, dt_da, DAtimestepsreq,  &
                       Obs )


!****************************************************************************************
!*                                                                                      *
!*  Output the observations structure                                                   *
!*                                                                                      *
!*  Filename          - Output filename                                                 *
!*  maxtime           - maximum time (seconds)                                          *
!*  dt                - model timestep (for information)                                *
!*  timestepsreq      - number of model timesteps required to cover all obs             *
!*  dt_da             - data assimilation timestep                                      *
!*  DAtimestepsreq    - number of data assimilation timesteps required to cover all obs *
!*  Obs               - Pointer to the first observation (linked list)                  *
!*                                                                                      *
!*  R. Bannister, 1.4da 08-04-2018                                                      *
!*                                                                                      *
!****************************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  Obs_type

IMPLICIT NONE

! Parameters
!-----------
CHARACTER(LEN=320),      INTENT(IN) :: Filename
INTEGER,                 INTENT(IN) :: maxtime
REAL(ZREAL8),            INTENT(IN) :: dt
INTEGER,                 INTENT(IN) :: timestepsreq
REAL(ZREAL8),            INTENT(IN) :: dt_da
INTEGER,                 INTENT(IN) :: DAtimestepsreq
TYPE(Obs_type), POINTER, INTENT(IN) :: Obs

! Local Variables
!----------------
TYPE(Obs_type), POINTER             :: thisob


  OPEN (13, file=Filename)
  WRITE (13,'(A)') 'Observation file for ABC model'
  WRITE (13,'(A18,I3)')    'Format version  : ', 1
  WRITE (13,'(A18,I8)')    'maxtime (s)     : ', maxtime
  WRITE (13,'(A18,E12.3)') 'Model ts (s)    : ', dt
  WRITE (13,'(A18,I8)')    'No model ts     : ', timestepsreq
  WRITE (13,'(A18,E12.3)') 'DA ts (s)       : ', dt_da
  WRITE (13,'(A18,I8)')    'No DA ts        : ', DAtimestepsreq
  WRITE (13,'(A)') '-----------------------------------------'

  thisob  => Obs

  DO
    IF (ASSOCIATED(thisob)) THEN
      IF (thisob % ob_ok) THEN
        WRITE (13,'(A)') '-----------------------------------------'
        WRITE (13,'(A18,I10)')   'Observation No  : ', thisob % obnumber_thisfile
        WRITE (13,'(A18,I6)')    'Batch ID        : ', thisob % batch
        WRITE (13,'(A18,I8)')    'Time of obs (s) : ', thisob % t
        WRITE (13,'(A18,F12.3)') 'Longitude (deg) : ', thisob % longitude_deg
        WRITE (13,'(A18,F12.3)') 'Height (m)      : ', thisob % level_ht
        WRITE (13,'(A18,I10)')   'xbox_lower      : ', thisob % xbox_lower
        WRITE (13,'(A18,I10)')   'xbox_lower_ws   : ', thisob % xbox_lower_ws
        WRITE (13,'(A18,I10)')   'zbox_lower      : ', thisob % zbox_lower
        WRITE (13,'(A18,I10)')   'zbox_lower_ws   : ', thisob % zbox_lower_ws
        WRITE (13,'(A18,I10)')   'tstep_lower     : ', thisob % tstep_lower


        WRITE (13,'(A18,I3)')    'Observation of  : ', thisob % ob_of_what
        SELECT CASE (thisob % ob_of_what)
        CASE (1)
          WRITE (13,'(A19)')    '                : u'
        CASE (2)
          WRITE (13,'(A19)')    '                : v'
        CASE (3)
          WRITE (13,'(A19)')    '                : w'
        CASE (4)
          WRITE (13,'(A19)')    '                : r'
        CASE (5)
          WRITE (13,'(A19)')    '                : b'
        CASE (6)
          WRITE (13,'(A24)')    '                : tracer'
        CASE (7)
          WRITE (13,'(A35)')    '                : horiz. wind speed'
        CASE (8)
          WRITE (13,'(A34)')    '                : total wind speed'
        END SELECT
        WRITE (13,'(A18,L)')     'y_true_known    : ', thisob % y_true_known
        WRITE (13,'(A18,E14.5)') 'y_true          : ', thisob % y_true
        WRITE (13,'(A18,E14.5)') 'y               : ', thisob % y
        WRITE (13,'(A18,E14.5)') 'stddev          : ', thisob % stddev
        WRITE (13,'(A18,E14.5)') 'y_ref           : ', thisob % y_ref
        WRITE (13,'(A18,E14.5)') 'd               : ', thisob % d
        WRITE (13,'(A18,E14.5)') 'deltay_m        : ', thisob % deltay_m
        WRITE (13,'(A18,E14.5)') 'hxmy            : ', thisob % hxmy
        WRITE (13,'(A18,E14.5)') 'deltay_m_hat    : ', thisob % deltay_m_hat
      END IF
      thisob => thisob % next
    ELSE
      EXIT
    END IF

  END DO

  CLOSE (13)

END SUBROUTINE Write_Obs
