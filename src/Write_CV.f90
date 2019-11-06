SUBROUTINE Write_CV (filename, CV_data, type_of_cv, CVT, longs, levs)

!********************************************************
!* Subroutine to write a control variable to a file     *
!*                                                      *
!*                                                      *
!* R. Bannister, vn1.4da, 11-12-17                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
  CV_type,                 &
  CVT_type,                &
  nlongs,                  &
  nlevs,                   &
  A, B, C, f


IMPLICIT NONE

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE '/usr/include/netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)     :: filename
TYPE(CV_type),    INTENT(IN)     :: CV_data
INTEGER,          INTENT(IN)     :: type_of_cv ! 1 = parameters
                                               ! 2 = stddev removed
                                               ! 3 = between spatial transforms
                                               ! 4 = actual control
                                               ! Add 10 for adjoint of above
TYPE(CVT_type),   INTENT(IN)     :: CVT        ! For options
REAL(ZREAL8),     INTENT(IN)     :: longs(1:nlongs)
REAL(ZREAL8),     INTENT(IN)     :: levs(1:nlevs)

!Declare local variables
!------------------------
INTEGER                      :: dimidScalar, varidScalar
INTEGER                      :: dimidLongs, varidLongs
INTEGER                      :: dimidLevs, varidLevs
INTEGER                      :: dimidwn, varidwn
INTEGER                      :: dimidvertmode, varidvertmode

INTEGER                      :: varid_A, varid_B, varid_C, varid_f, varid_typeofcv
INTEGER                      :: varid_CVT_order, varid_CVT_param_opt_gb, varid_CVT_param_opt_hb
INTEGER                      :: varid_CVT_param_opt_ab, varid_CVT_param_opt_reg
INTEGER                      :: varid_CVT_vert_opt_sym, varid_CVT_stddev_opt
INTEGER                      :: varid_v1, varid_v2, varid_v3, varid_v4, varid_v5, varid_v6

INTEGER                      :: ncid, ierr, ddA(1), ddB(2)
INTEGER                      :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7, ierr8, ierr9, ierr10, ierr11, ierr12
INTEGER                      :: startA(1:1), countA(1:1)
INTEGER                      :: startB(1:2), countB(1:2)

REAL(ZREAL8)                 :: wavenumbers(1:nlongs)
REAL(ZREAL8)                 :: vert_mode(1:nlevs)
INTEGER                      :: scalar(1:1)
INTEGER                      :: x, l, k
INTEGER                      :: real_index, imag_index

!*****************************************************************************************
!PRINT*, 'Write control variable'
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
  ierr1 = NF_DEF_DIM(ncid, 'scalar', 1, dimidScalar)
  ierr2 = NF_DEF_DIM(ncid, 'longs', nlongs, dimidLongs)
  ierr3 = NF_DEF_DIM(ncid, 'level', nlevs, dimidLevs)
  ierr4 = NF_DEF_DIM(ncid, 'wavenumber', nlongs, dimidwn)
  ierr5 = NF_DEF_DIM(ncid, 'vert_mode', nlevs, dimidvertmode)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5) /= 0) THEN
     PRINT*, '***Error defining dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids defined'
  END IF

  !Define the variables (include variables giving the dim. values)
  !---------------------------------------------------------------

  ! Dimension variables

  ddA(1) = dimidScalar
  ierr1  = NF_DEF_VAR(ncid, 'scalar', NF_INT, 1, ddA, varidScalar)

  ddA(1) = dimidLongs
  ierr2  = NF_DEF_VAR(ncid, 'longs', NF_DOUBLE, 1, ddA, varidLongs)

  ddA(1) = dimidLevs
  ierr3  = NF_DEF_VAR(ncid, 'level', NF_DOUBLE, 1, ddA, varidLevs)

  ddA(1) = dimidwn
  ierr4  = NF_DEF_VAR(ncid, 'wavenumber', NF_DOUBLE, 1, ddA, varidwn)

  ddA(1) = dimidvertmode
  ierr5  = NF_DEF_VAR(ncid, 'vert_mode', NF_DOUBLE, 1, ddA, varidvertmode)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr5) /= 0) THEN
     PRINT*, '***Error defining dimension variable ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     STOP
  !ELSE
    !PRINT*, 'Dimension variable ids defined'
  END IF

  ! Main variables

  !--- Parameters and options ---
  ddA(1)  = dimidScalar
  ierr1   = NF_DEF_VAR(ncid, 'A', NF_DOUBLE, 1, ddA, varid_A)
  ierr2   = NF_DEF_VAR(ncid, 'B', NF_DOUBLE, 1, ddA, varid_B)
  ierr3   = NF_DEF_VAR(ncid, 'C', NF_DOUBLE, 1, ddA, varid_C)
  ierr4   = NF_DEF_VAR(ncid, 'f', NF_DOUBLE, 1, ddA, varid_f)
  ierr5   = NF_DEF_VAR(ncid, 'type_of_cv', NF_INT, 1, ddA, varid_typeofcv)
  ierr6   = NF_DEF_VAR(ncid, 'cvt_order', NF_INT, 1, ddA, varid_CVT_order)
  ierr7   = NF_DEF_VAR(ncid, 'cvt_param_opt_gb', NF_INT, 1, ddA, varid_CVT_param_opt_gb)
  ierr8   = NF_DEF_VAR(ncid, 'cvt_param_opt_hb', NF_INT, 1, ddA, varid_CVT_param_opt_hb)
  ierr9   = NF_DEF_VAR(ncid, 'cvt_param_opt_ab', NF_INT, 1, ddA, varid_CVT_param_opt_ab)
  ierr10  = NF_DEF_VAR(ncid, 'cvt_param_opt_reg', NF_INT, 1, ddA, varid_CVT_param_opt_reg)
  ierr11  = NF_DEF_VAR(ncid, 'cvt_vert_opt_sym', NF_INT, 1, ddA, varid_CVT_vert_opt_sym)
  ierr12  = NF_DEF_VAR(ncid, 'cvt_stddev_opt', NF_INT, 1, ddA, varid_CVT_stddev_opt)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11 + ierr12) /= 0) THEN
    PRINT*, '***Error defining main variable ids (batch 1)***'
    PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
    PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
    PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
    PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
    PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
    PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
    PRINT*, 'ierr7 ',  ierr7,  NF_STRERROR(ierr7)
    PRINT*, 'ierr8 ',  ierr8,  NF_STRERROR(ierr8)
    PRINT*, 'ierr9 ',  ierr9,  NF_STRERROR(ierr9)
    PRINT*, 'ierr10',  ierr10, NF_STRERROR(ierr10)
    PRINT*, 'ierr11',  ierr11, NF_STRERROR(ierr11)
    PRINT*, 'ierr12',  ierr12, NF_STRERROR(ierr12)
    STOP
  END IF

  ! --- Main fields --
  SELECT CASE (CVT % CVT_order)
  CASE (1)
    ! As original MetO order

    SELECT CASE (CVT % CVT_vert_opt_sym)
    CASE (1)
      ! Original MetO order, non-symmetric vertical transform

      SELECT CASE (type_of_cv)
      CASE (1,11) ! parameters
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (2,12) ! stddev removed
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (3,13) ! between spatial transforms
        ddB(1) = dimidLongs
        ddB(2) = dimidvertmode
      CASE (4,14) ! actual control
        ddB(1) = dimidwn
        ddB(2) = dimidvertmode
      END SELECT

    CASE (2)
      ! Original MetO order, symmetric vertical transform

      SELECT CASE (type_of_cv)
      CASE (1,11) ! parameters
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (2,12) ! stddev removed
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (3,13) ! between spatial transforms
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (4,14) ! actual control
        ddB(1) = dimidwn
        ddB(2) = dimidLevs
      END SELECT

    END SELECT

  CASE (2)
    ! Reversed MetO order

    SELECT CASE (CVT % CVT_vert_opt_sym)
    CASE (1)
      ! Reversed MetO order, non-symmetric vertical transform

      SELECT CASE (type_of_cv)
      CASE (1,11) ! parameters
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (2,12) ! stddev removed
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (3,13) ! between spatial transforms
        ddB(1) = dimidwn
        ddB(2) = dimidLevs
      CASE (4,14) ! actual control
        ddB(1) = dimidwn
        ddB(2) = dimidvertmode
      END SELECT

    CASE (2)
      ! Original MetO order, symmetric vertical transform

      SELECT CASE (type_of_cv)
      CASE (1,11) ! parameters
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (2,12) ! stddev removed
        ddB(1) = dimidLongs
        ddB(2) = dimidLevs
      CASE (3,13) ! between spatial transforms
        ddB(1) = dimidwn
        ddB(2) = dimidLevs
      CASE (4,14) ! actual control
        ddB(1) = dimidwn
        ddB(2) = dimidLevs
      END SELECT

    END SELECT

  CASE (3)
    ! As REP thesis
    PRINT *, 'Normal mode-based transform not implemented'
    STOP

  END SELECT

  ierr1   = NF_DEF_VAR(ncid, 'v_1', NF_DOUBLE, 2, ddB, varid_v1)
  ierr2   = NF_DEF_VAR(ncid, 'v_2', NF_DOUBLE, 2, ddB, varid_v2)
  ierr3   = NF_DEF_VAR(ncid, 'v_3', NF_DOUBLE, 2, ddB, varid_v3)
  ierr4   = NF_DEF_VAR(ncid, 'v_4', NF_DOUBLE, 2, ddB, varid_v4)
  ierr5   = NF_DEF_VAR(ncid, 'v_5', NF_DOUBLE, 2, ddB, varid_v5)
  ierr6   = NF_DEF_VAR(ncid, 'v_6', NF_DOUBLE, 2, ddB, varid_v6)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
    PRINT*, '***Error defining main variable ids (batch 2)***'
    PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
    PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
    PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
    PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
    PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
    PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
    STOP
  END IF



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



  ! --------------------------------------------
  ! Output the values of the dimension variables
  ! --------------------------------------------
  scalar(1) = 0
  startA(1) = 1
  countA(1) = 1
  ierr1     = NF_PUT_VARA_INT(ncid, varidScalar, startA, countA, scalar(1:1))
  countA(1) = nlongs
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varidLongs, startA, countA, longs(1:nlongs))

  countA(1) = nlevs
  ierr3     = NF_PUT_VARA_DOUBLE(ncid, varidLevs, startA, countA, levs(1:nlevs))

  wavenumbers(1) = 0.0
  DO k = 2, nlongs/2
    real_index = 2*k-2
    imag_index = 2*k-1
    wavenumbers(real_index) = REAL(k-1)
    wavenumbers(imag_index) = REAL(k-1) + 0.5  !Actually the imag part
  END DO
  wavenumbers(nlongs) = nlongs/2
  countA(1) = nlongs
  ierr4     = NF_PUT_VARA_DOUBLE(ncid, varidwn, startA, countA, wavenumbers(1:nlongs))

  DO l = 1, nlevs
    vert_mode(l) = REAL(l)
  END DO
  countA(1) = nlevs
  ierr5     = NF_PUT_VARA_DOUBLE(ncid, varidvertmode, startA, countA, vert_mode(1:nlevs))

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr5) /= 0) THEN
    PRINT*, '***Error outputting dimension values ***'
    PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
    PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
    PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
    PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
    PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
    STOP
  END IF


  !--------------------------------------------
  ! Output the values of the main variables
  ! -------------------------------------------

  !--- Parameters and options ---

  ierr1     = NF_PUT_VAR1_DOUBLE(ncid, varid_A, 1, A)
  ierr2     = NF_PUT_VAR1_DOUBLE(ncid, varid_B, 1, B)
  ierr3     = NF_PUT_VAR1_DOUBLE(ncid, varid_C, 1, C)
  ierr4     = NF_PUT_VAR1_DOUBLE(ncid, varid_f, 1, f)
  ierr5     = NF_PUT_VAR1_INT(ncid, varid_typeofcv, 1, type_of_cv)
  ierr6     = NF_PUT_VAR1_INT(ncid, varid_CVT_order, 1, CVT % CVT_order)
  ierr7     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_gb, 1, CVT % CVT_param_opt_gb)
  ierr8     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_hb, 1, CVT % CVT_param_opt_hb)
  ierr9     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_ab, 1, CVT % CVT_param_opt_ab)
  ierr10    = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_reg, 1, CVT % CVT_param_opt_reg)
  ierr11    = NF_PUT_VAR1_INT(ncid, varid_CVT_vert_opt_sym, 1, CVT % CVT_vert_opt_sym)
  ierr12    = NF_PUT_VAR1_INT(ncid, varid_CVT_stddev_opt, 1, CVT % CVT_stddev_opt)
  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11 + ierr12) /= 0) THEN
    PRINT*, '***Error writing main variables (batch 1)***'
    PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
    PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
    PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
    PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
    PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
    PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
    PRINT*, 'ierr7 ',  ierr7,  NF_STRERROR(ierr7)
    PRINT*, 'ierr8 ',  ierr8,  NF_STRERROR(ierr8)
    PRINT*, 'ierr9 ',  ierr9,  NF_STRERROR(ierr9)
    PRINT*, 'ierr10',  ierr10, NF_STRERROR(ierr10)
    PRINT*, 'ierr11',  ierr11, NF_STRERROR(ierr11)
    PRINT*, 'ierr12',  ierr12, NF_STRERROR(ierr12)
    STOP
  END IF


  ! --- Control fields themselves ---
  startB(1) = 1
  countB(1) = nlongs
  startB(2) = 1
  countB(2) = nlevs
  ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_v1, startB, countB, CV_data % v1(1:nlongs, 1:nlevs))
  ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_v2, startB, countB, CV_data % v2(1:nlongs, 1:nlevs))
  ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_v3, startB, countB, CV_data % v3(1:nlongs, 1:nlevs))
  ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_v4, startB, countB, CV_data % v4(1:nlongs, 1:nlevs))
  ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_v5, startB, countB, CV_data % v5(1:nlongs, 1:nlevs))
  ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_v6, startB, countB, CV_data % v6(1:nlongs, 1:nlevs))

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
    PRINT*, '***Error writing main variables (batch 2)***'
    PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
    PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
    PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
    PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
    PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
    PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
    STOP
  END IF


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


END SUBROUTINE Write_CV
