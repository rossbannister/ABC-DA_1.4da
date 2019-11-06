SUBROUTINE Write_one_field (filename, Nx, Nz, field, fieldname)

!********************************************************
!* Subroutine to write a single 2d field                *
!*                                                      *
!*                                                      *
!* R. Bannister, vn1.4da, 14-12-17                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8

IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE '/usr/include/netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)  :: filename
INTEGER,          INTENT(IN)  :: Nx, Nz
REAL(ZREAL8),     INTENT(IN)  :: field(1:Nx, 1:Nz)
CHARACTER(LEN=*), INTENT(IN)  :: fieldname


!Declare local variables
!------------------------
INTEGER                       :: ncid, ierr, ddA(1), ddB(2)
INTEGER                       :: dimidx, dimidz
INTEGER                       :: varidx, varidz, varid_field
INTEGER                       :: startA(1), countA(1), startB(2), countB(2)
INTEGER                       :: ierr1, ierr2
REAL(ZREAL8)                  :: xs(1:Nx)
REAL(ZREAL8)                  :: zs(1:Nz)
INTEGER                       :: x, z

!*****************************************************************************************
!PRINT*, 'Write one field'
!*****************************************************************************************

  ! Create netCDF file
  !-------------------------------------
  ierr = NF_CREATE(filename, NF_CLOBBER, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error creating file ***'
    PRINT*, filename
    PRINT*, ierr, NF_STRERROR(ierr)
    STOP
    !ELSE
    !PRINT*, 'FILE CREATED'
  ENDIF

  !Define the dimensions
  !---------------------
  ierr1 = NF_DEF_DIM(ncid, 'x', Nx, dimidx)
  ierr2 = NF_DEF_DIM(ncid, 'z', Nz, dimidz)

  IF ((ierr1 + ierr2) /= 0) THEN
     PRINT*, '***Error defining dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids defined'
  ENDIF

  !Define the variables (include variables giving the dim. values)
  !---------------------------------------------------------------

  ! Dimension variables

  ddA(1) = dimidx
  ierr1  = NF_DEF_VAR(ncid, 'x', NF_DOUBLE, 1, ddA, varidx)

  ddA(1) = dimidz
  ierr2  = NF_DEF_VAR(ncid, 'z', NF_DOUBLE, 1, ddA, varidz)

  IF ((ierr1 + ierr2) /= 0) THEN
     PRINT*, '***Error defining dimension variable ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     STOP
  !ELSE
  !PRINT*, 'Dimension variable ids defined'
  ENDIF

  ! Main variable

  ddB(1)  = dimidx
  ddB(2)  = dimidz
  ierr    = NF_DEF_VAR(ncid, fieldname, NF_DOUBLE, 2, ddB, varid_field)

  IF (ierr /= 0) THEN
    PRINT*, '***Error defining main variable id ***'
    PRINT*, 'ierr ',  ierr,  NF_STRERROR(ierr)
    STOP
  ENDIF

  !Change mode of netCDF operation from define to write
  !------------------------------------------------------
  ierr = NF_ENDDEF(ncid)

  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error changing mode of netCDF operation ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  !ELSE
  !PRINT*, 'Mode Changed'
  END IF

  !---------------------------------------------
  ! Output the values of the dimension variables
  ! --------------------------------------------
  DO x = 1, Nx
    xs(x) = REAL(x)
  END DO
  DO z = 1, Nz
    zs(z) = REAL(z)
  END DO

  startA(1) = 1
  countA(1) = Nx
  ierr1     = NF_PUT_VARA_DOUBLE(ncid, varidx, startA, countA, xs(1:Nx))
  countA(1) = Nz
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varidz, startA, countA, zs(1:Nz))

  IF ((ierr1 + ierr2) /= 0) THEN
    PRINT*, '***Error writing dimension data ***'
    PRINT*, 'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*, 'ierr2', ierr2, NF_STRERROR(ierr2)
    STOP
  !ELSE
    !PRINT*, 'Dimension Variables output ok'
  ENDIF

!--------------------------------------------
! Output the values of the main variable
! -------------------------------------------

startB(1) = 1
countB(1) = Nx
startB(2) = 1
countB(2) = Nz

ierr      = NF_PUT_VARA_DOUBLE(ncid, varid_field, startB, countB, field(1:Nx, 1:Nz))

IF (ierr /= 0) THEN
  PRINT*, '***Error writing main variable data ***'
  PRINT*,'ierr',   ierr,   NF_STRERROR(ierr)
  STOP
!ELSE
  !PRINT*, 'Main data written'
ENDIF

!Close-up the file
!-----------------
ierr = NF_CLOSE(ncid)

IF ( ierr .NE. 0 ) THEN
  PRINT*, ' *** Error closing netCDF file ***'
  PRINT*,'ierr', ierr, NF_STRERROR(ierr)
  PRINT*,'xconv -i ', filename, ' &'
  STOP
!ELSE
  !PRINT*, 'File closed'
ENDIF


END SUBROUTINE Write_one_field
