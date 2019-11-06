SUBROUTINE Read_state_2d (filename, state, dims, t_input)

!********************************************************
!* Subroutine to read a field of 2d toy model data      *
!*        of ABC_type                                   *
!*                                                      *
!* t is the read time (-1 for last time in the file)    *
!*                                                      *
!* R. Petrie, 2.0 10-6-11                               *
!* R. Petrie, 3.0 30-7-13                               *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
  ABC_type,                &
  dims_type,               &
  nlongs,                  &
  nlevs

IMPLICIT NONE

INCLUDE "Boundaries.interface"

!NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE '/usr/include/netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*),   INTENT(IN)    :: filename
TYPE(ABC_type),     INTENT(INOUT) :: state
TYPE(dims_type),    INTENT(INOUT) :: dims
INTEGER,            INTENT(IN)    :: t_input

!Declare local variables
!-----------------------
INTEGER         :: ncid, ierr, dd(3)
INTEGER         :: dimidLongs_u, varidLongs_u
INTEGER         :: dimidLongs_v, varidLongs_v
INTEGER         :: dimidHalf_levs, varidHalf_levs
INTEGER         :: dimidFull_levs, varidFull_levs
INTEGER         :: dimid_time, varid_time
INTEGER         :: dimidEnergy, varidEnergy
INTEGER         :: varid_u, varid_v, varid_w
INTEGER         :: varid_b, varid_beff
INTEGER         :: varid_r, varid_rho
INTEGER         :: varid_tracer, varid_hydro, varid_geost
INTEGER         :: varid_ke, varid_be, varid_ee, varid_te
INTEGER         :: startA(1), countA(1), startC(3), countC(3)
INTEGER         :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7
INTEGER         :: ierr8, ierr9, ierr10, ierr11, ierr12, ierr13, ierr14
INTEGER         :: ierr15, ierr16
INTEGER         :: z, ntimes, t
REAL(ZREAL8)    :: temp(nlongs), tempvert(nlevs)

!PRINT*,'read field called'

!Open the netCDF file
!---------------------
ierr = NF_OPEN(filename, NF_NOWRITE, ncid)
IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error opening file ***'
  PRINT*, ierr, NF_STRERROR(ierr)
  PRINT*, 'FILE :: ', filename
  STOP
ENDIF

!Get the dimension ids
!---------------------
ierr1 = NF_INQ_DIMID(ncid, 'longs_u', dimidLongs_u)
ierr2 = NF_INQ_DIMID(ncid, 'longs_v', dimidLongs_v)
ierr3 = NF_INQ_DIMID(ncid, 'half_level', dimidHalf_levs)
ierr4 = NF_INQ_DIMID(ncid, 'full_level', dimidFull_levs)
ierr5 = NF_INQ_DIMID(ncid, 'time', dimid_time)

IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5) /= 0) THEN
  PRINT*, '***Error getting dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  STOP
!ELSE
!  PRINT*, ' Dimension ids ok'
ENDIF

! Find out how many times are in the file
ierr = NF_INQ_DIMLEN (ncid, dimid_time, ntimes)
IF (ierr /= 0) THEN
  PRINT*, '***Error getting length of time dimension ***'
  PRINT*,'ierr', ierr, NF_STRERROR(ierr)
  STOP
END IF

! Check or set the time of interest
IF ((t_input > ntimes) .OR. (t_input == 0) .OR. (t_input < -1)) THEN
  ! Out of range
  PRINT *, 'Error reading in time: ', t_input, ' is out of range'
  STOP
END IF
IF (t_input == -1) THEN
  ! Set the time of interest to the last one in the file
  t = ntimes
ELSE
  ! Set the time of interest to that specified
  t = t_input
END IF
PRINT *, 'Retrieving data for time ', t

! Get the variable ids, for dimensions
!-------------------------------------
ierr1 = NF_INQ_VARID(ncid, 'longs_u', varidLongs_u)
ierr2 = NF_INQ_VARID(ncid, 'longs_v', varidLongs_v)
ierr3 = NF_INQ_VARID(ncid, 'half_level', varidHalf_levs)
ierr4 = NF_INQ_VARID(ncid, 'full_level', varidFull_levs)

IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
  PRINT*, '***Error getting variable dimension ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
!  PRINT*, ' Variable dimension ids ok'
ENDIF


!Get the main variable ids
!-------------------------
ierr1  = NF_INQ_VARID(ncid, 'u', varid_u)
ierr2  = NF_INQ_VARID(ncid, 'v', varid_v)
ierr3  = NF_INQ_VARID(ncid, 'w', varid_w)
ierr4  = NF_INQ_VARID(ncid, 'r_prime', varid_r)
ierr5  = NF_INQ_VARID(ncid, 'b_prime', varid_b)
ierr6  = NF_INQ_VARID(ncid, 'rho', varid_rho)
ierr7  = NF_INQ_VARID(ncid, 'b_effective', varid_beff)
ierr8  = NF_INQ_VARID(ncid, 'tracer', varid_tracer)
ierr9  = NF_INQ_VARID(ncid, 'geo_imbal', varid_geost)
ierr10 = NF_INQ_VARID(ncid, 'hydro_imbal', varid_hydro)
ierr11 = NF_INQ_VARID(ncid, 'ke', varid_ke)
ierr12 = NF_INQ_VARID(ncid, 'be', varid_be)
ierr13 = NF_INQ_VARID(ncid, 'ee', varid_ee)
ierr14 = NF_INQ_VARID(ncid, 'te', varid_te)

IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
      ierr9 + ierr10 + ierr11 + ierr12 + ierr13 + ierr14 ) /= 0) THEN
  PRINT*, '***Error getting variable ids ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
  PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
  PRINT*,'ierr7', ierr7, NF_STRERROR(ierr7)
  PRINT*,'ierr8', ierr8, NF_STRERROR(ierr8)
  PRINT*,'ierr9', ierr9, NF_STRERROR(ierr9)
  PRINT*,'ierr10', ierr10, NF_STRERROR(ierr10)
  PRINT*,'ierr11', ierr11, NF_STRERROR(ierr11)
  PRINT*,'ierr12', ierr12, NF_STRERROR(ierr12)
  PRINT*,'ierr13', ierr13, NF_STRERROR(ierr13)
  PRINT*,'ierr14', ierr14, NF_STRERROR(ierr14)
  STOP
ELSE
!  PRINT*, ' Variable ids are ok'
ENDIF


!Get the dimension data from the file
!------------------------------------
startA(1) = 1
countA(1) = nlongs
ierr1 = NF_GET_VARA_DOUBLE (ncid, varidLongs_u, startA, countA, dims % longs_u(1:nlongs))
ierr2 = NF_GET_VARA_DOUBLE (ncid, varidLongs_v, startA, countA, dims % longs_v(1:nlongs))

countA(1) = nlevs
ierr3 = NF_GET_VARA_DOUBLE (ncid, varidHalf_levs, startA, countA, dims % half_levs(1:nlevs))
ierr4 = NF_GET_VARA_DOUBLE (ncid, varidFull_levs, startA, countA, dims % full_levs(1:nlevs))

dims % half_levs(0) = - dims % half_levs(1)
dims % full_levs(0) = 0.0
dims % half_levs(nlevs+1) = 2. * dims % half_levs(nlevs) - dims % half_levs(nlevs - 1)
dims % full_levs(nlevs+1) = 2. * dims % full_levs(nlevs) - dims % full_levs(nlevs - 1)

IF ((ierr1 + ierr2 + ierr3 + ierr4) /= 0) THEN
  PRINT*, '*** Error getting dimension data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
!ELSE
!  PRINT*, ' Dimension data ok'
ENDIF


!Get the main data from the file
!-------------------------------
startC(1) = 1
countC(1) = nlongs !longs
startC(2) = 1
countC(2) = 1      !levs
startC(3) = t
countC(3) = 1      !time

! Read 2D variables
DO z = 1, nlevs

  startC(2) = z

  ! Read u
  !-------
  ierr1 = NF_GET_VARA_DOUBLE (ncid, varid_u, startC, countC,  temp)
  state % u(1:nlongs,z) = temp(1:nlongs)

  ! Read v
  !-------
  ierr2 = NF_GET_VARA_DOUBLE (ncid, varid_v, startC, countC, temp)
  state % v(1:nlongs, z) = temp(1:nlongs)

  ! Read w
  !-------
  ierr3 = NF_GET_VARA_DOUBLE (ncid, varid_w, startC, countC, temp)
  state % w(1:nlongs, z) = temp(1:nlongs)

  ! Read r
  !--------
  ierr4 = NF_GET_VARA_DOUBLE (ncid, varid_r, startC, countC, temp)
  state % r(1:nlongs, z) = temp(1:nlongs)

  ! Read bp
  !--------
  ierr5 = NF_GET_VARA_DOUBLE (ncid, varid_b, startC, countC, temp)
  state % b(1:nlongs, z) = temp(1:nlongs)

  ! Read rho
  !---------
  ierr6 = NF_GET_VARA_DOUBLE (ncid, varid_rho, startC, countC, temp)
  state % rho (1:nlongs,z) = temp(1:nlongs)

  ! Read beffective
  !----------------
  ierr7 = NF_GET_VARA_DOUBLE (ncid, varid_beff, startC, countC, temp)
  state % b_ef (1:nlongs,z) = temp(1:nlongs)

  ! Read tracer
  !------------
  ierr8 = NF_GET_VARA_DOUBLE (ncid, varid_tracer, startC, countC, temp)
  state % tracer(1:nlongs, z) = temp(1:nlongs)

  ! Read geostrohpic imbalance
  !---------------------------
  ierr9 = NF_GET_VARA_DOUBLE (ncid, varid_geost, startC, countC, temp)
  state % geost_imbal(1:nlongs, z) = temp(1:nlongs)

  ! Read hydrostatic imbalance
  !---------------------------
  ierr10 = NF_GET_VARA_DOUBLE (ncid, varid_hydro, startC, countC, temp)
  state % hydro_imbal(1:nlongs, z) = temp(1:nlongs)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + &
       ierr9 + ierr10  ) /= 0) THEN
    PRINT*, '*** Error getting main variable data ***'
    PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
    PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
    PRINT*,'ierr7', ierr7, NF_STRERROR(ierr7)
    PRINT*,'ierr8', ierr8, NF_STRERROR(ierr8)
    PRINT*,'ierr9', ierr9, NF_STRERROR(ierr9)
    STOP
  ENDIF

ENDDO

! Read energy
!------------
startA(1) = t
countA(1) = 1
ierr1 = NF_GET_VARA_DOUBLE (ncid, varid_ke, startA, countA, state % Kinetic_Energy)
ierr2 = NF_GET_VARA_DOUBLE (ncid, varid_be, startA, countA, state % Buoyant_Energy)
ierr3 = NF_GET_VARA_DOUBLE (ncid, varid_ee, startA, countA, state % Elastic_Energy)
ierr4 = NF_GET_VARA_DOUBLE (ncid, varid_te, startA, countA, state % Total_Energy)

IF (( ierr1 + ierr2 + ierr3 + ierr4 ) /= 0) THEN
  PRINT*, '*** Error getting main variable data ***'
  PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
  PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
  PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
  PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
  STOP
ENDIF

!Close the netCDF file
ierr = NF_CLOSE(ncid)
IF ( ierr .NE. 0 ) THEN
   PRINT*, '***ERROR closing file***'
   STOP
ENDIF

CALL Boundaries (state)

END SUBROUTINE Read_state_2d
