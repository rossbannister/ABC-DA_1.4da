SUBROUTINE Read_CV (filename, CV_data, type_of_cv, CVT, &
                    read_only_opts, consistency_check)

!********************************************************
!* Subroutine to read a control variable from a file    *
!*                                                      *
!*                                                      *
!* R. Bannister, vn1.4da, 14-12-17                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
  CV_type,                 &
  CVT_type,                &
  nlongs,                  &
  nlevs,                   &
  A, B, C, f,              &
  small


IMPLICIT NONE

INCLUDE "Boundaries_CV.interface"

! NetCDF library (file format used to read/write data)
!----------------------------------------------------
INCLUDE '/usr/include/netcdf.inc'

!Declare parameters
!------------------
CHARACTER(LEN=*), INTENT(IN)     :: filename
TYPE(CV_type),    INTENT(INOUT)  :: CV_data
INTEGER,          INTENT(IN)     :: type_of_cv ! 1 = parameters
                                               ! 2 = stddev removed
                                               ! 3 = between spatial transforms
                                               ! 4 = actual control
                                               ! Add 10 for adjoint of above
TYPE(CVT_type),   INTENT(IN)     :: CVT        ! For options
LOGICAL,          INTENT(IN)     :: read_only_opts       ! Read only A, etc, and opts
LOGICAL,          INTENT(IN)     :: consistency_check    ! Do consistency check

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

INTEGER                      :: ncid, ierr
INTEGER                      :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7, ierr8, ierr9, ierr10, ierr11, ierr12
INTEGER                      :: startA(1:1), countA(1:1)
INTEGER                      :: startB(1:2), countB(1:2)

INTEGER                      :: type_of_cv_in
REAL(ZREAL8)                 :: A_in, B_in, C_in, f_in


!*****************************************************************************************
!PRINT*, 'Read control variable'
!*****************************************************************************************

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
  ierr1 = NF_INQ_DIMID(ncid, 'scalar', dimidScalar)
  IF (.NOT.read_only_opts) THEN
    ierr2 = NF_INQ_DIMID(ncid, 'longs', dimidLongs)
    ierr3 = NF_INQ_DIMID(ncid, 'level', dimidLevs)
    ierr4 = NF_INQ_DIMID(ncid, 'wavenumber', dimidwn)
    ierr5 = NF_INQ_DIMID(ncid, 'vert_mode', dimidvertmode)
  ELSE
    ierr2 = 0
    ierr3 = 0
    ierr4 = 0
    ierr5 = 0
  END IF

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5) /= 0) THEN
     PRINT*, '***Error getting dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     STOP
  !ELSE
  !   PRINT*, 'Dimension ids ok'
  END IF

  ! Get the variable ids, for dimensions
  !-------------------------------------
  ierr1  = NF_INQ_VARID(ncid, 'scalar', varidScalar)
  IF (.NOT.read_only_opts) THEN
    ierr2  = NF_INQ_VARID(ncid, 'longs', varidLongs)
    ierr4  = NF_INQ_VARID(ncid, 'level', varidLevs)
    ierr5  = NF_INQ_VARID(ncid, 'wavenumber', varidwn)
    ierr6  = NF_INQ_VARID(ncid, 'vert_mode', varidvertmode)
  ELSE
    ierr2 = 0
    ierr3 = 0
    ierr4 = 0
    ierr5 = 0
  END IF

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr5 + ierr6) /= 0) THEN
     PRINT*, '***Error getting variable dimension ids ***'
     PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
     PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
     PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
     PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
     PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
     PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
     STOP
  !ELSE
    !PRINT*, 'Dimension variable ids ok'
  END IF

  !Get the main variable ids
  !-------------------------

  !--- Parameters and options ---
  ierr1   = NF_INQ_VARID(ncid, 'A', varid_A)
  ierr2   = NF_INQ_VARID(ncid, 'B', varid_B)
  ierr3   = NF_INQ_VARID(ncid, 'C', varid_C)
  ierr4   = NF_INQ_VARID(ncid, 'f', varid_f)
  ierr5   = NF_INQ_VARID(ncid, 'type_of_cv', varid_typeofcv)
  ierr6   = NF_INQ_VARID(ncid, 'cvt_order', varid_CVT_order)
  ierr7   = NF_INQ_VARID(ncid, 'cvt_param_opt_gb', varid_CVT_param_opt_gb)
  ierr8   = NF_INQ_VARID(ncid, 'cvt_param_opt_hb', varid_CVT_param_opt_hb)
  ierr9   = NF_INQ_VARID(ncid, 'cvt_param_opt_ab', varid_CVT_param_opt_ab)
  ierr10  = NF_INQ_VARID(ncid, 'cvt_param_opt_reg', varid_CVT_param_opt_reg)
  ierr11  = NF_INQ_VARID(ncid, 'cvt_vert_opt_sym', varid_CVT_vert_opt_sym)
  ierr12  = NF_INQ_VARID(ncid, 'cvt_stddev_opt', varid_CVT_stddev_opt)

  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11 + ierr12) /= 0) THEN
    PRINT*, '***Error getting main variable ids (batch 1)***'
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


  IF (.NOT.read_only_opts) THEN
    ! --- Main fields --
    ierr1   = NF_INQ_VARID(ncid, 'v_1', varid_v1)
    ierr2   = NF_INQ_VARID(ncid, 'v_2', varid_v2)
    ierr3   = NF_INQ_VARID(ncid, 'v_3', varid_v3)
    ierr4   = NF_INQ_VARID(ncid, 'v_4', varid_v4)
    ierr5   = NF_INQ_VARID(ncid, 'v_5', varid_v5)
    ierr6   = NF_INQ_VARID(ncid, 'v_6', varid_v6)

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error getting main variable ids (batch 2)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF
  END IF


  ! --------------------------------------------
  ! Input the values of the parameters and options
  ! --------------------------------------------
  ierr1     = NF_GET_VAR1_DOUBLE(ncid, varid_A, 1, A_in)
  ierr2     = NF_GET_VAR1_DOUBLE(ncid, varid_B, 1, B_in)
  ierr3     = NF_GET_VAR1_DOUBLE(ncid, varid_C, 1, C_in)
  ierr4     = NF_GET_VAR1_DOUBLE(ncid, varid_f, 1, f_in)
  ierr5     = NF_GET_VAR1_INT(ncid, varid_typeofcv, 1, type_of_cv_in)
  ierr6     = NF_GET_VAR1_INT(ncid, varid_CVT_order, 1, CVT % CVT_order)
  ierr7     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_gb, 1, CVT % CVT_param_opt_gb)
  ierr8     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_hb, 1, CVT % CVT_param_opt_hb)
  ierr9     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_ab, 1, CVT % CVT_param_opt_ab)
  ierr10    = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_reg, 1, CVT % CVT_param_opt_reg)
  ierr11    = NF_GET_VAR1_INT(ncid, varid_CVT_vert_opt_sym, 1, CVT % CVT_vert_opt_sym)
  ierr12    = NF_GET_VAR1_INT(ncid, varid_CVT_stddev_opt, 1, CVT % CVT_stddev_opt)
  IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11 + ierr12) /= 0) THEN
    PRINT*, '***Error getting main variables (batch 1)***'
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


  IF (.NOT.read_only_opts) THEN
    ! --------------------------------------------
    ! Input the values of the main variables
    ! --------------------------------------------

    startB(1) = 1
    countB(1) = nlongs
    startB(2) = 1
    countB(2) = nlevs
    ierr1 = NF_GET_VARA_DOUBLE (ncid, varid_v1, startB, countB,   &
                                CV_data % v1(1:nlongs, 1:nlevs))
    ierr2 = NF_GET_VARA_DOUBLE (ncid, varid_v2, startB, countB,   &
                                CV_data % v2(1:nlongs, 1:nlevs))
    ierr3 = NF_GET_VARA_DOUBLE (ncid, varid_v3, startB, countB,   &
                                CV_data % v3(1:nlongs, 1:nlevs))
    ierr4 = NF_GET_VARA_DOUBLE (ncid, varid_v4, startB, countB,   &
                                CV_data % v4(1:nlongs, 1:nlevs))
    ierr5 = NF_GET_VARA_DOUBLE (ncid, varid_v5, startB, countB,   &
                                CV_data % v5(1:nlongs, 1:nlevs))
    ierr6 = NF_GET_VARA_DOUBLE (ncid, varid_v6, startB, countB,   &
                                CV_data % v6(1:nlongs, 1:nlevs))

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error getting main variables (batch 1)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

  END IF

  !Close-up the file
  !-----------------
  ierr = NF_CLOSE(ncid)

  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error closing netCDF file ***'
    PRINT*,'ierr', ierr, NF_STRERROR(ierr)
    STOP
  !ELSE
    !PRINT*, 'File closed'
  ENDIF


  ! Do a consistency check of the parameters
  IF (consistency_check .AND. (.NOT.read_only_opts)) THEN
    IF ((type_of_cv_in /= type_of_cv) .OR. &
        (DABS(A_in - A) > small)      .OR. &
        (DABS(B_in - B) > small)      .OR. &
        (DABS(C_in - C) > small)      .OR. &
        (DABS(f_in - f) > small)) THEN
      PRINT *, 'Mismatch(es) detected in input control variable file'
      PRINT *, filename
      PRINT *, 'In file value (type_of_cv):', type_of_cv_in
      PRINT *, 'Expecting                 :', type_of_cv
      PRINT *, 'In file value (A)         :', A_in
      PRINT *, 'Expecting                 :', A
      PRINT *, 'In file value (B)         :', B_in
      PRINT *, 'Expecting                 :', B
      PRINT *, 'In file value (C)         :', C_in
      PRINT *, 'Expecting                 :', C
      PRINT *, 'In file value (f)         :', f_in
      PRINT *, 'Expecting                 :', f
      STOP
    END IF
  END IF


END SUBROUTINE Read_CV
