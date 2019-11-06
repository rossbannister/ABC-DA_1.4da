SUBROUTINE Read_um_data_2d (umdata, filename, latitude)

!**********************************
!* Subroutine to read a latitude  *
!* slice of um netcdf output into *
!* a compound data type um_data   *
!*                                *
!* R. Petrie                      *
!* version 2                      *
!* 06/06/2011                     *
!**********************************

USE DefConsTypes, ONLY :  &
    ZREAL8,               &
    UM_type,              &
    nlongs, nlevs

IMPLICIT NONE

INCLUDE '/usr/include/netcdf.inc'

!Declare parameters
!------------------
TYPE(UM_type),        INTENT(INOUT) :: umdata
CHARACTER (LEN=*),    INTENT(IN)    :: filename
INTEGER,              INTENT(IN)    :: latitude

! Declare local variables
!-------------------------
INTEGER             :: ncid, ierr
INTEGER             :: dimidLongs_u, dimidLongs_v, dimidhalf_levs, dimidfull_levs
INTEGER             :: varidLongs_u, varidLongs_v, varidhalf_levs, varidfull_levs
INTEGER             :: varidu, varidv, varidw, variddensity, varidtheta
INTEGER             :: varidorog, varidexpres
INTEGER             :: startA(1), countA(1), startB(4), countB(4),z
INTEGER             :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
REAL (ZREAL8)       :: temp(360)

! Open the netCDF file
!----------------------
ierr = NF_OPEN(TRIM(filename), NF_NOWRITE, ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, 'Filename: ', TRIM(filename)
  PRINT*, ierr, NF_STRERROR(ierr)
  STOP
ENDIF

!Get the dimension ids
!-----------------------
ierr = NF_INQ_DIMID(ncid, 'x', dimidLongs_u)
ierr1 = ierr
ierr = NF_INQ_DIMID(ncid, 'x_1', dimidLongs_v)
ierr2 = ierr
ierr = NF_INQ_DIMID(ncid, 'hybrid_ht_3', dimidhalf_levs)
ierr3 = ierr
ierr = NF_INQ_DIMID(ncid, 'hybrid_ht_2', dimidfull_levs)
ierr4 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
! PRINT*, ' Dimension ids ok'
ENDIF

! Get the dimension variable ids
!---------------------------------
ierr = NF_INQ_VARID(ncid, 'x', varidLongs_u)
ierr1 = ierr
ierr = NF_INQ_VARID(ncid, 'x_1', varidLongs_v)
ierr2 = ierr
ierr = NF_INQ_VARID(ncid, 'hybrid_ht_3', varidhalf_levs)
ierr3 = ierr
ierr = NF_INQ_VARID(ncid, 'hybrid_ht_2', varidfull_levs)
ierr4 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
!ELSE
!  PRINT*, ' Dimension variable ids ok'
ENDIF


! Get the variable ids
!----------------------
ierr = NF_INQ_VARID(ncid, 'u', varidu)
ierr1 = ierr
ierr = NF_INQ_VARID(ncid, 'v', varidv)
ierr2 = ierr
ierr = NF_INQ_VARID(ncid, 'dz_dt', varidw)
ierr3 = ierr
ierr = NF_INQ_VARID(ncid, 'unspecified', variddensity)
ierr4 = ierr
ierr = NF_INQ_VARID(ncid, 'theta', varidtheta)
ierr5 = ierr
ierr = NF_INQ_VARID(ncid, 'ht', varidorog)
ierr6 = ierr
ierr = NF_INQ_VARID(ncid, 'field7', varidexpres)
ierr7 = ierr
IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. (ierr4 .NE. 0) .OR. (ierr5 .NE. 0)&
      .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
  PRINT*, '***Error getting variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
  STOP
!ELSE
!  PRINT*, ' Variable ids ok'
ENDIF

! Get the dimension data from the file
!---------------------------------------
! Longitudes
!-------------
startA(1) = 1
countA(1) = nlongs

ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, umdata % longs_u(1:nlongs))
ierr2 = ierr
ierr = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, umdata % longs_v(1:nlongs))
ierr3 = ierr

! Level Heights
!----------------
startA(1) = 1
countA(1) = nlevs+1
ierr = NF_GET_VARA_DOUBLE (ncid, varidfull_levs, startA, countA, umdata % full_levs(0:nlevs))
ierr4 = ierr
ierr = NF_GET_VARA_DOUBLE (ncid, varidhalf_levs, startA, countA, umdata % half_levs(1:nlevs+1))
ierr5 = ierr

IF ( (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
     (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) ) THEN
  PRINT*, '***Error getting dimension data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
! ELSE
!  PRINT*, ' Dimension data ok'
ENDIF

! Get the main data from the file
!----------------------------------
startB(1) = 1
countB(1) = 360       ! All longitude points
startB(2) = latitude  ! Selected latitude slice
countB(2) = 1         !lats
startB(3) = 1
countB(3) = 1         !levs
startB(4) = 1
countB(4) = 1         !time

! u
!----
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidu, startB, countB,  temp)
  umdata % u(1:nlongs,z) = temp(1:nlongs)
END DO
ierr1 = ierr

! v
!----
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidv, startB, countB, temp)
  umdata % v(1:nlongs, z) = temp(1:nlongs)
END DO
ierr2 = ierr

! w
!----
DO z=1, nlevs+1
  startB(3) = z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidw, startB, countB, temp)
  umdata % w(1:nlongs, z-1) = temp(1:nlongs)
ENDDO
ierr3 = ierr

! density
!---------
DO z=1, nlevs
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, variddensity, startB, countB, temp)
  umdata % density(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr4 = ierr

! theta
!-------
DO z=1, nlevs
  startB(3)= z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidtheta, startB, countB, temp)
  umdata % theta(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr5 = ierr

! exner presure
!----------------
DO z=1, nlevs + 1
  startB(3) = z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidexpres, startB, countB, temp)
  umdata % exner_pressure(1:nlongs, z) = temp(1:nlongs)
ENDDO
ierr6 = ierr

! orographic height
!--------------------
DO z=1,1
  startB(3)=z
  ierr = NF_GET_VARA_DOUBLE (ncid, varidorog, startB, countB, umdata % orog_height(1:nlongs))
ENDDO
ierr7 = ierr

IF ( (ierr1 .NE. 0) .OR. (ierr2 .NE. 0) .OR. (ierr3 .NE. 0) .OR. &
     (ierr4 .NE. 0) .OR. (ierr5 .NE. 0) .OR. (ierr6 .NE. 0) .OR. (ierr7 .NE. 0) ) THEN
  PRINT*, '***Error getting main variable data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr6, NF_STRERROR(ierr6)
  STOP
! ELSE
!  PRINT*, ' Main variable data ok'
ENDIF

!Close the netCDF file
ierr = NF_CLOSE(ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, '***ERROR closing file***'
  STOP
ENDIF

END SUBROUTINE Read_um_data_2d
