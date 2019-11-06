SUBROUTINE Write_Covs (filename, CVT, longs, levs,  &
                       newfile, output_params,      &        ! flags
                       output_stddevs, output_vert, &        ! flags
                       output_horiz, output_regression)      ! flags

!********************************************************
!* Subroutine to write an ABC cov structure to a file   *
!*                                                      *
!* R. Bannister, vn1.4da, 08-12-17                      *
!*                                                      *
!********************************************************

USE DefConsTypes, ONLY :   &
  ZREAL8,                  &
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
CHARACTER(LEN=*),  INTENT(IN) :: filename
TYPE(CVT_type),    INTENT(IN) :: CVT
REAL(ZREAL8),      INTENT(IN) :: longs(1:nlongs)
REAL(ZREAL8),      INTENT(IN) :: levs(1:nlevs)
LOGICAL, OPTIONAL, INTENT(IN) :: newfile
LOGICAL, OPTIONAL, INTENT(IN) :: output_params
LOGICAL, OPTIONAL, INTENT(IN) :: output_stddevs
LOGICAL, OPTIONAL, INTENT(IN) :: output_vert
LOGICAL, OPTIONAL, INTENT(IN) :: output_horiz
LOGICAL, OPTIONAL, INTENT(IN) :: output_regression


!Declare local variables
!------------------------
INTEGER                       :: dimidScalar, varidScalar
INTEGER                       :: dimidLongs, varidLongs
INTEGER                       :: dimidLevs, varidLevs
INTEGER                       :: dimidwn, varidwn
INTEGER                       :: dimidvertmode, varidvertmode

LOGICAL                       :: set_newfile, set_output_params, set_output_stddevs, set_output_vert
LOGICAL                       :: set_output_horiz, set_output_regression

INTEGER                       :: varid_A, varid_B, varid_C, varid_f
INTEGER                       :: varid_CVT_order, varid_CVT_param_opt_gb, varid_CVT_param_opt_hb
INTEGER                       :: varid_CVT_param_opt_ab, varid_CVT_param_opt_reg
INTEGER                       :: varid_CVT_vert_opt_sym, varid_CVT_stddev_opt
INTEGER                       :: varid_sigma1, varid_sigma2, varid_sigma3, varid_sigma4, varid_sigma5, varid_sigma6
INTEGER                       :: varid_VertMode1, varid_VertMode2, varid_VertMode3, varid_VertMode4
INTEGER                       :: varid_VertMode5, varid_VertMode6
INTEGER                       :: varid_VertEV1, varid_VertEV2, varid_VertEV3, varid_VertEV4, varid_VertEV5, varid_VertEV6
INTEGER                       :: varid_HorizEV1, varid_HorizEV2, varid_HorizEV3, varid_HorizEV4, varid_HorizEV5, varid_HorizEV6
INTEGER                       :: varid_Regression, varid_Cov_rbalrbal, varid_Cov_rtotrbal

INTEGER                       :: ncid, ierr, i, ddA(1), ddB(2), ddC(3)
INTEGER                       :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7, ierr8, ierr9, ierr10, ierr11
INTEGER                       :: startA(1:1), countA(1:1)
INTEGER                       :: startB(1:2), countB(1:2)
INTEGER                       :: startC(1:3), countC(1:3)

REAL(ZREAL8)                  :: wavenumbers(1:nlongs/2+1)
REAL(ZREAL8)                  :: vert_mode(1:nlevs)
INTEGER                       :: scalar(1:1)
INTEGER                       :: x, l, k



!*****************************************************************************************
!PRINT*, 'Write Covs'
!*****************************************************************************************

  IF (PRESENT(newfile)) THEN
    set_newfile = newfile
  ELSE
    set_newfile = .FALSE.
  END IF

  IF (PRESENT(output_params)) THEN
    set_output_params = output_params
  ELSE
    set_output_params = .FALSE.
  END IF

  IF (PRESENT(output_stddevs)) THEN
    set_output_stddevs = output_stddevs
  ELSE
    set_output_stddevs = .FALSE.
  END IF

  IF (PRESENT(output_vert)) THEN
    set_output_vert = output_vert
  ELSE
    set_output_vert = .FALSE.
  END IF

  IF (PRESENT(output_horiz)) THEN
    set_output_horiz = output_horiz
  ELSE
    set_output_horiz = .FALSE.
  END IF

  IF (PRESENT(output_regression)) THEN
    set_output_regression = output_regression
  ELSE
    set_output_regression = .FALSE.
  END IF

  IF (set_newfile) THEN
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
    END IF


    !Define the dimensions
    !---------------------
    ierr1 = NF_DEF_DIM(ncid, 'scalar', 1, dimidScalar)
    ierr2 = NF_DEF_DIM(ncid, 'longs', nlongs, dimidLongs)
    ierr3 = NF_DEF_DIM(ncid, 'level', nlevs, dimidLevs)
    ierr4 = NF_DEF_DIM(ncid, 'wavenumber', nlongs/2+1, dimidwn)
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
    ierr4  = NF_DEF_VAR(ncid, 'level', NF_DOUBLE, 1, ddA, varidLevs)

    ddA(1) = dimidwn
    ierr5  = NF_DEF_VAR(ncid, 'wavenumber', NF_DOUBLE, 1, ddA, varidwn)

    ddA(1) = dimidvertmode
    ierr6  = NF_DEF_VAR(ncid, 'vert_mode', NF_DOUBLE, 1, ddA, varidvertmode)

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr5 + ierr6) /= 0) THEN
       PRINT*, '***Error defining dimension variable ids ***'
       PRINT*,'ierr1', ierr1, NF_STRERROR(ierr1)
       PRINT*,'ierr2', ierr2, NF_STRERROR(ierr2)
       PRINT*,'ierr3', ierr3, NF_STRERROR(ierr3)
       PRINT*,'ierr4', ierr4, NF_STRERROR(ierr4)
       PRINT*,'ierr5', ierr5, NF_STRERROR(ierr5)
       PRINT*,'ierr6', ierr6, NF_STRERROR(ierr6)
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
    ierr5   = NF_DEF_VAR(ncid, 'cvt_order', NF_INT, 1, ddA, varid_CVT_order)
    ierr6   = NF_DEF_VAR(ncid, 'cvt_param_opt_gb', NF_INT, 1, ddA, varid_CVT_param_opt_gb)
    ierr7   = NF_DEF_VAR(ncid, 'cvt_param_opt_hb', NF_INT, 1, ddA, varid_CVT_param_opt_hb)
    ierr8   = NF_DEF_VAR(ncid, 'cvt_param_opt_ab', NF_INT, 1, ddA, varid_CVT_param_opt_ab)
    ierr9   = NF_DEF_VAR(ncid, 'cvt_param_opt_reg', NF_INT, 1, ddA, varid_CVT_param_opt_reg)
    ierr10  = NF_DEF_VAR(ncid, 'cvt_vert_opt_sym', NF_INT, 1, ddA, varid_CVT_vert_opt_sym)
    ierr11  = NF_DEF_VAR(ncid, 'cvt_stddev_opt', NF_INT, 1, ddA, varid_CVT_stddev_opt)

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11) /= 0) THEN
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
      STOP
    END IF

    ! --- Standard deviations --
    SELECT CASE (CVT % CVT_stddev_opt)
    CASE (1)
      ! Stddev constant for each control variable
      ddA(1) = dimidScalar
      ierr1   = NF_DEF_VAR(ncid, 'sigma1', NF_DOUBLE, 1, ddA, varid_sigma1)
      ierr2   = NF_DEF_VAR(ncid, 'sigma2', NF_DOUBLE, 1, ddA, varid_sigma2)
      ierr3   = NF_DEF_VAR(ncid, 'sigma3', NF_DOUBLE, 1, ddA, varid_sigma3)
      ierr4   = NF_DEF_VAR(ncid, 'sigma4', NF_DOUBLE, 1, ddA, varid_sigma4)
      ierr5   = NF_DEF_VAR(ncid, 'sigma5', NF_DOUBLE, 1, ddA, varid_sigma5)
      ierr6   = NF_DEF_VAR(ncid, 'sigma6', NF_DOUBLE, 1, ddA, varid_sigma6)

    CASE (2)
      ! Stddev level dependent
      ddA(1)  = dimidLevs
      ierr1   = NF_DEF_VAR(ncid, 'sigma1', NF_DOUBLE, 1, ddA, varid_sigma1)
      ierr2   = NF_DEF_VAR(ncid, 'sigma2', NF_DOUBLE, 1, ddA, varid_sigma2)
      ierr3   = NF_DEF_VAR(ncid, 'sigma3', NF_DOUBLE, 1, ddA, varid_sigma3)
      ierr4   = NF_DEF_VAR(ncid, 'sigma4', NF_DOUBLE, 1, ddA, varid_sigma4)
      ierr5   = NF_DEF_VAR(ncid, 'sigma5', NF_DOUBLE, 1, ddA, varid_sigma5)
      ierr6   = NF_DEF_VAR(ncid, 'sigma6', NF_DOUBLE, 1, ddA, varid_sigma6)

    CASE (3)
      ! Stddev longitude and level dependent
      ddB(1)  = dimidLongs
      ddB(2)  = dimidLevs
      ierr1   = NF_DEF_VAR(ncid, 'sigma1', NF_DOUBLE, 2, ddB, varid_sigma1)
      ierr2   = NF_DEF_VAR(ncid, 'sigma2', NF_DOUBLE, 2, ddB, varid_sigma2)
      ierr3   = NF_DEF_VAR(ncid, 'sigma3', NF_DOUBLE, 2, ddB, varid_sigma3)
      ierr4   = NF_DEF_VAR(ncid, 'sigma4', NF_DOUBLE, 2, ddB, varid_sigma4)
      ierr5   = NF_DEF_VAR(ncid, 'sigma5', NF_DOUBLE, 2, ddB, varid_sigma5)
      ierr6   = NF_DEF_VAR(ncid, 'sigma6', NF_DOUBLE, 2, ddB, varid_sigma6)

    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error defining main variable ids (batch 2)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF


    ! --- Vertical modes ---
    SELECT CASE (CVT % CVT_order)

    CASE (1)
      ! Original MetO order of transforms (vert modes constant with longitude)
      ddB(1)  = dimidLevs
      ddB(2)  = dimidLevs
      ierr1   = NF_DEF_VAR(ncid, 'vertmode1', NF_DOUBLE, 2, ddB, varid_VertMode1)
      ierr2   = NF_DEF_VAR(ncid, 'vertmode2', NF_DOUBLE, 2, ddB, varid_VertMode2)
      ierr3   = NF_DEF_VAR(ncid, 'vertmode3', NF_DOUBLE, 2, ddB, varid_VertMode3)
      ierr4   = NF_DEF_VAR(ncid, 'vertmode4', NF_DOUBLE, 2, ddB, varid_VertMode4)
      ierr5   = NF_DEF_VAR(ncid, 'vertmode5', NF_DOUBLE, 2, ddB, varid_VertMode5)
      ierr6   = NF_DEF_VAR(ncid, 'vertmode6', NF_DOUBLE, 2, ddB, varid_VertMode6)

    CASE (2)
      ! Reversed MetO order of transforms (vert modes fn of horiz wavenumber)
      ddC(1)  = dimidLevs
      ddC(2)  = dimidLevs
      ddC(3)  = dimidwn
      ierr1   = NF_DEF_VAR(ncid, 'vertmode1', NF_DOUBLE, 3, ddC, varid_VertMode1)
      ierr2   = NF_DEF_VAR(ncid, 'vertmode2', NF_DOUBLE, 3, ddC, varid_VertMode2)
      ierr3   = NF_DEF_VAR(ncid, 'vertmode3', NF_DOUBLE, 3, ddC, varid_VertMode3)
      ierr4   = NF_DEF_VAR(ncid, 'vertmode4', NF_DOUBLE, 3, ddC, varid_VertMode4)
      ierr5   = NF_DEF_VAR(ncid, 'vertmode5', NF_DOUBLE, 3, ddC, varid_VertMode5)
      ierr6   = NF_DEF_VAR(ncid, 'vertmode6', NF_DOUBLE, 3, ddC, varid_VertMode6)
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error defining main variable ids (batch 3)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF


    ! --- Vertical eigenvalues ---
    SELECT CASE (CVT % CVT_order)

    CASE (1)
      ! Original MetO order of transforms (eigenvalues constant with longitude)
      ddA(1)  = dimidvertmode
      ierr1   = NF_DEF_VAR(ncid, 'vertEV1', NF_DOUBLE, 1, ddA, varid_VertEV1)
      ierr2   = NF_DEF_VAR(ncid, 'vertEV2', NF_DOUBLE, 1, ddA, varid_VertEV2)
      ierr3   = NF_DEF_VAR(ncid, 'vertEV3', NF_DOUBLE, 1, ddA, varid_VertEV3)
      ierr4   = NF_DEF_VAR(ncid, 'vertEV4', NF_DOUBLE, 1, ddA, varid_VertEV4)
      ierr5   = NF_DEF_VAR(ncid, 'vertEV5', NF_DOUBLE, 1, ddA, varid_VertEV5)
      ierr6   = NF_DEF_VAR(ncid, 'vertEV6', NF_DOUBLE, 1, ddA, varid_VertEV6)

    CASE (2)
      ! Reversed MetO order of transforms (eigenvalues fn of horiz wavenumber)
      ddB(1)  = dimidvertmode
      ddB(2)  = dimidwn
      ierr1   = NF_DEF_VAR(ncid, 'vertEV1', NF_DOUBLE, 2, ddB, varid_VertEV1)
      ierr2   = NF_DEF_VAR(ncid, 'vertEV2', NF_DOUBLE, 2, ddB, varid_VertEV2)
      ierr3   = NF_DEF_VAR(ncid, 'vertEV3', NF_DOUBLE, 2, ddB, varid_VertEV3)
      ierr4   = NF_DEF_VAR(ncid, 'vertEV4', NF_DOUBLE, 2, ddB, varid_VertEV4)
      ierr5   = NF_DEF_VAR(ncid, 'vertEV5', NF_DOUBLE, 2, ddB, varid_VertEV5)
      ierr6   = NF_DEF_VAR(ncid, 'vertEV6', NF_DOUBLE, 2, ddB, varid_VertEV6)
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error defining main variable ids (batch 4)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF


    ! --- Horizontal eigenvalues ---
    IF ((CVT % CVT_order == 1) .AND. (CVT % CVT_vert_opt_sym == 1)) THEN
      ! Original MetO order and non-symmetric vert trans (wn and vert mode)
      ddB(1)  = dimidwn
      ddB(2)  = dimidvertmode
    ELSE
      ! Original MetO order and symmetric vert trans OR
      ! Reversed MetO order (wn and vert level)
      ddB(1)  = dimidwn
      ddB(2)  = dimidLevs
    END IF
    ierr1   = NF_DEF_VAR(ncid, 'horizEV1', NF_DOUBLE, 2, ddB, varid_HorizEV1)
    ierr2   = NF_DEF_VAR(ncid, 'horizEV2', NF_DOUBLE, 2, ddB, varid_HorizEV2)
    ierr3   = NF_DEF_VAR(ncid, 'horizEV3', NF_DOUBLE, 2, ddB, varid_HorizEV3)
    ierr4   = NF_DEF_VAR(ncid, 'horizEV4', NF_DOUBLE, 2, ddB, varid_HorizEV4)
    ierr5   = NF_DEF_VAR(ncid, 'horizEV5', NF_DOUBLE, 2, ddB, varid_HorizEV5)
    ierr6   = NF_DEF_VAR(ncid, 'horizEV6', NF_DOUBLE, 2, ddB, varid_HorizEV6)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error defining main variable ids (batch 5)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF


    ! --- Vertical regression matrix for balanced pressure (and other vert covs) ---
    IF (CVT % CVT_param_opt_reg == 1) THEN
      ddB(1)  = dimidLevs
      ddB(2)  = dimidLevs
      ierr1   = NF_DEF_VAR(ncid, 'regress_gb', NF_DOUBLE, 2, ddB, varid_Regression)
      ierr2   = NF_DEF_VAR(ncid, 'cov_rbalrbal', NF_DOUBLE, 2, ddB, varid_Cov_rbalrbal)
      ierr3   = NF_DEF_VAR(ncid, 'cov_rtotrbal', NF_DOUBLE, 2, ddB, varid_Cov_rtotrbal)
      IF ((ierr1 + ierr2 + ierr3 ) /= 0) THEN
        PRINT*, '***Error defining main variable ids (batch 6)***'
        PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
        PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
        PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
        STOP
      END IF
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

    DO k = 1, nlongs/2+1
      wavenumbers(k) = REAL(k-1)
    END DO
    countA(1) = nlongs/2+1
    ierr4     = NF_PUT_VARA_DOUBLE(ncid, varidwn, startA, countA, wavenumbers(1:nlongs/2+1))

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




  ELSE

    ! Open existing netCDF file for output
    !-------------------------------------
    ierr = NF_OPEN(filename, NF_WRITE, ncid)
    IF ( ierr .NE. 0 ) THEN
      PRINT*, ' *** Error reopening file ***'
      PRINT*, filename
      PRINT*, ierr, NF_STRERROR(ierr)
      STOP
      !ELSE
      !PRINT*, 'FILE REOPENED'
    END IF



  END IF


  !-----------------------------------------------------
  ! Output the values of the main variables as requested
  ! ----------------------------------------------------

  IF (set_output_params) THEN


    !--- Parameters and options ---
    ierr1     = NF_INQ_VARID(ncid, 'A', varid_A)
    ierr2     = NF_INQ_VARID(ncid, 'B', varid_B)
    ierr3     = NF_INQ_VARID(ncid, 'C', varid_C)
    ierr4     = NF_INQ_VARID(ncid, 'f', varid_f)
    ierr5     = NF_INQ_VARID(ncid, 'cvt_order', varid_CVT_order)
    ierr6     = NF_INQ_VARID(ncid, 'cvt_param_opt_gb', varid_CVT_param_opt_gb)
    ierr7     = NF_INQ_VARID(ncid, 'cvt_param_opt_hb', varid_CVT_param_opt_hb)
    ierr8     = NF_INQ_VARID(ncid, 'cvt_param_opt_ab', varid_CVT_param_opt_ab)
    ierr9     = NF_INQ_VARID(ncid, 'cvt_param_opt_reg', varid_CVT_param_opt_reg)
    ierr10    = NF_INQ_VARID(ncid, 'cvt_vert_opt_sym', varid_CVT_vert_opt_sym)
    ierr11    = NF_INQ_VARID(ncid, 'cvt_stddev_opt', varid_CVT_stddev_opt)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11) /= 0) THEN
      PRINT*, '***Error equiring about variables (batch 1)***'
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
      STOP
    END IF
    startA(1) = 1
    countA(1) = 1
    ierr1     = NF_PUT_VAR1_DOUBLE(ncid, varid_A, 1, A)
    ierr2     = NF_PUT_VAR1_DOUBLE(ncid, varid_B, 1, B)
    ierr3     = NF_PUT_VAR1_DOUBLE(ncid, varid_C, 1, C)
    ierr4     = NF_PUT_VAR1_DOUBLE(ncid, varid_f, 1, f)
    ierr5     = NF_PUT_VAR1_INT(ncid, varid_CVT_order, 1, CVT % CVT_order)
    ierr6     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_gb, 1, CVT % CVT_param_opt_gb)
    ierr7     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_hb, 1, CVT % CVT_param_opt_hb)
    ierr8     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_ab, 1, CVT % CVT_param_opt_ab)
    ierr9     = NF_PUT_VAR1_INT(ncid, varid_CVT_param_opt_reg, 1, CVT % CVT_param_opt_reg)
    ierr10    = NF_PUT_VAR1_INT(ncid, varid_CVT_vert_opt_sym, 1, CVT % CVT_vert_opt_sym)
    ierr11    = NF_PUT_VAR1_INT(ncid, varid_CVT_stddev_opt, 1, CVT % CVT_stddev_opt)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11) /= 0) THEN
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
      STOP
    END IF

  END IF



  IF (set_output_stddevs) THEN

    ! --- Standard deviations --

    ierr1     = NF_INQ_VARID(ncid, 'sigma1', varid_sigma1)
    ierr2     = NF_INQ_VARID(ncid, 'sigma2', varid_sigma2)
    ierr3     = NF_INQ_VARID(ncid, 'sigma3', varid_sigma3)
    ierr4     = NF_INQ_VARID(ncid, 'sigma4', varid_sigma4)
    ierr5     = NF_INQ_VARID(ncid, 'sigma5', varid_sigma5)
    ierr6     = NF_INQ_VARID(ncid, 'sigma6', varid_sigma6)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error equiring about variables (batch 2)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

    SELECT CASE (CVT % CVT_stddev_opt)
    CASE (1)
      ! Stddev constant for each control variable
      ierr1   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma1, 1, CVT % sigma1(1,1))
      ierr2   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma2, 1, CVT % sigma2(1,1))
      ierr3   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma3, 1, CVT % sigma3(1,1))
      ierr4   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma4, 1, CVT % sigma4(1,1))
      ierr5   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma5, 1, CVT % sigma5(1,1))
      ierr6   = NF_PUT_VAR1_DOUBLE(ncid, varid_sigma6, 1, CVT % sigma6(1,1))

    CASE (2)
      ! Stddev level dependent
      startA(1) = 1
      countA(1) = nlevs
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma1, startA, countA, CVT % sigma1(1, 1:nlevs))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma2, startA, countA, CVT % sigma2(1, 1:nlevs))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma3, startA, countA, CVT % sigma3(1, 1:nlevs))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma4, startA, countA, CVT % sigma4(1, 1:nlevs))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma5, startA, countA, CVT % sigma5(1, 1:nlevs))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma6, startA, countA, CVT % sigma6(1, 1:nlevs))

    CASE (3)
      ! Stddev longitude and level dependent
      startB(1) = 1
      countB(1) = nlongs
      startB(2) = 1
      countB(2) = nlevs
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma1, startB, countB, CVT % sigma1(1:nlongs, 1:nlevs))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma2, startB, countB, CVT % sigma2(1:nlongs, 1:nlevs))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma3, startB, countB, CVT % sigma3(1:nlongs, 1:nlevs))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma4, startB, countB, CVT % sigma4(1:nlongs, 1:nlevs))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma5, startB, countB, CVT % sigma5(1:nlongs, 1:nlevs))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_sigma6, startB, countB, CVT % sigma6(1:nlongs, 1:nlevs))

    END SELECT

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

  END IF


  IF (set_output_vert) THEN

    ! --- Vertical modes ---

    ierr1     = NF_INQ_VARID(ncid, 'vertmode1', varid_VertMode1)
    ierr2     = NF_INQ_VARID(ncid, 'vertmode2', varid_VertMode2)
    ierr3     = NF_INQ_VARID(ncid, 'vertmode3', varid_VertMode3)
    ierr4     = NF_INQ_VARID(ncid, 'vertmode4', varid_VertMode4)
    ierr5     = NF_INQ_VARID(ncid, 'vertmode5', varid_VertMode5)
    ierr6     = NF_INQ_VARID(ncid, 'vertmode6', varid_VertMode6)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error equiring about variables (batch 3)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

    SELECT CASE (CVT % CVT_order)

    CASE (1)
      ! Original MetO order of transforms (vert modes constant with longitude)
      startB(1) = 1
      countB(1) = nlevs
      startB(2) = 1
      countB(2) = nlevs
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode1, startB, countB, CVT % VertMode1(1:nlevs, 1:nlevs, 1))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode2, startB, countB, CVT % VertMode2(1:nlevs, 1:nlevs, 1))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode3, startB, countB, CVT % VertMode3(1:nlevs, 1:nlevs, 1))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode4, startB, countB, CVT % VertMode4(1:nlevs, 1:nlevs, 1))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode5, startB, countB, CVT % VertMode5(1:nlevs, 1:nlevs, 1))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode6, startB, countB, CVT % VertMode6(1:nlevs, 1:nlevs, 1))

    CASE (2)
      ! Reversed MetO order of transforms (vert modes fn of horiz wavenumber)
      startC(1) = 1
      countC(1) = nlevs
      startC(2) = 1
      countC(2) = nlevs
      startC(3) = 1
      countC(3) = nlongs/2+1
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode1, startC, countC, CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode2, startC, countC, CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode3, startC, countC, CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode4, startC, countC, CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode5, startC, countC, CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_VertMode6, startC, countC, CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1))
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error writing main variables (batch 3)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF


    ! --- Vertical eigenvalues ---
    ierr1     = NF_INQ_VARID(ncid, 'vertEV1', varid_VertEV1)
    ierr2     = NF_INQ_VARID(ncid, 'vertEV2', varid_VertEV2)
    ierr3     = NF_INQ_VARID(ncid, 'vertEV3', varid_VertEV3)
    ierr4     = NF_INQ_VARID(ncid, 'vertEV4', varid_VertEV4)
    ierr5     = NF_INQ_VARID(ncid, 'vertEV5', varid_VertEV5)
    ierr6     = NF_INQ_VARID(ncid, 'vertEV6', varid_VertEV6)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error equiring about variables (batch 4)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF
    SELECT CASE (CVT % CVT_order)

    CASE (1)
      ! Original MetO order of transforms (eigenvalues constant with longitude)
      startA(1) = 1
      countA(1) = nlevs
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV1, startA, countA, CVT % VertEV1(1:nlevs, 1))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV2, startA, countA, CVT % VertEV2(1:nlevs, 1))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV3, startA, countA, CVT % VertEV3(1:nlevs, 1))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV4, startA, countA, CVT % VertEV4(1:nlevs, 1))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV5, startA, countA, CVT % VertEV5(1:nlevs, 1))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV6, startA, countA, CVT % VertEV6(1:nlevs, 1))

    CASE (2)
      ! Reversed MetO order of transforms (eigenvalues fn of horiz wavenumber)
      startB(1) = 1
      countB(1) = nlevs
      startB(2) = 1
      countB(2) = nlongs/2+1
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV1, startB, countB, CVT % VertEV1(1:nlevs, 1:nlongs/2+1))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV2, startB, countB, CVT % VertEV2(1:nlevs, 1:nlongs/2+1))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV3, startB, countB, CVT % VertEV3(1:nlevs, 1:nlongs/2+1))
      ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV4, startB, countB, CVT % VertEV4(1:nlevs, 1:nlongs/2+1))
      ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV5, startB, countB, CVT % VertEV5(1:nlevs, 1:nlongs/2+1))
      ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_VertEV6, startB, countB, CVT % VertEV6(1:nlevs, 1:nlongs/2+1))
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error writing main variables (batch 4)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

  END IF


  IF (set_output_horiz) THEN
    ! --- Horizontal eigenvalues ---

    ierr1     = NF_INQ_VARID(ncid, 'horizEV1', varid_HorizEV1)
    ierr2     = NF_INQ_VARID(ncid, 'horizEV2', varid_HorizEV2)
    ierr3     = NF_INQ_VARID(ncid, 'horizEV3', varid_HorizEV3)
    ierr4     = NF_INQ_VARID(ncid, 'horizEV4', varid_HorizEV4)
    ierr5     = NF_INQ_VARID(ncid, 'horizEV5', varid_HorizEV5)
    ierr6     = NF_INQ_VARID(ncid, 'horizEV6', varid_HorizEV6)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error equiring about variables (batch 2)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

    startB(1) = 1
    countB(1) = nlongs/2+1
    startB(2) = 1
    countB(2) = nlevs
    ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV1, startB, countB, CVT % HorizEV1(1:nlongs/2+1, 1:nlevs))
    ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV2, startB, countB, CVT % HorizEV2(1:nlongs/2+1, 1:nlevs))
    ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV3, startB, countB, CVT % HorizEV3(1:nlongs/2+1, 1:nlevs))
    ierr4     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV4, startB, countB, CVT % HorizEV4(1:nlongs/2+1, 1:nlevs))
    ierr5     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV5, startB, countB, CVT % HorizEV5(1:nlongs/2+1, 1:nlevs))
    ierr6     = NF_PUT_VARA_DOUBLE(ncid, varid_HorizEV6, startB, countB, CVT % HorizEV6(1:nlongs/2+1, 1:nlevs))
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error writing main variables (batch 5)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF
  END IF


  IF (set_output_regression) THEN

    ! --- Vertical regression matrix for balanced pressure (and other vert covs) ---
    IF (CVT % CVT_param_opt_reg == 1) THEN
      ierr1     = NF_INQ_VARID(ncid, 'regress_gb', varid_Regression)
      ierr2     = NF_INQ_VARID(ncid, 'cov_rbalrbal', varid_Cov_rbalrbal)
      ierr3     = NF_INQ_VARID(ncid, 'cov_rtotrbal', varid_Cov_rtotrbal)
      IF ((ierr1 + ierr2 + ierr3 ) /= 0) THEN
        PRINT*, '***Error equiring about variables (batch 6)***'
        PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
        PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
        PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
        STOP
      END IF
      startB(1) = 1
      countB(1) = nlevs
      startB(2) = 1
      countB(2) = nlevs
      ierr1     = NF_PUT_VARA_DOUBLE(ncid, varid_Regression, startB, countB, CVT % Regression(1:nlevs, 1:nlevs))
      ierr2     = NF_PUT_VARA_DOUBLE(ncid, varid_Cov_rbalrbal, startB, countB, CVT % Cov_rbalrbal(1:nlevs, 1:nlevs))
      ierr3     = NF_PUT_VARA_DOUBLE(ncid, varid_Cov_rtotrbal, startB, countB, CVT % Cov_rtotrbal(1:nlevs, 1:nlevs))
      IF ((ierr1 + ierr2 + ierr3 ) /= 0) THEN
        PRINT*, '***Error writing main variables (batch 6)***'
        PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
        PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
        PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
        STOP
      END IF
    END IF
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
  END IF


END SUBROUTINE Write_Covs
