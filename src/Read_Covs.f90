SUBROUTINE Read_Covs (filename, CVT,                &
                      input_params,                 &     ! flags
                      input_stddevs, input_vert,    &     ! flags
                      input_horiz, input_regression)      ! flags

!********************************************************
!* Subroutine to read an ABC cov structure from a file  *
!*                                                      *
!* R. Bannister, vn1.4da, 16-12-17                      *
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
CHARACTER(LEN=*), INTENT(IN)    :: filename
TYPE(CVT_type),   INTENT(INOUT) :: CVT
LOGICAL,          INTENT(IN)    :: input_params
LOGICAL,          INTENT(IN)    :: input_stddevs
LOGICAL,          INTENT(IN)    :: input_vert
LOGICAL,          INTENT(IN)    :: input_horiz
LOGICAL,          INTENT(IN)    :: input_regression


!Declare local variables
!------------------------
INTEGER                      :: varid_A, varid_B, varid_C, varid_f
INTEGER                      :: varid_CVT_order, varid_CVT_param_opt_gb, varid_CVT_param_opt_hb
INTEGER                      :: varid_CVT_param_opt_ab, varid_CVT_param_opt_reg
INTEGER                      :: varid_CVT_vert_opt_sym, varid_CVT_stddev_opt
INTEGER                      :: varid_sigma1, varid_sigma2, varid_sigma3, varid_sigma4, varid_sigma5, varid_sigma6
INTEGER                      :: varid_VertMode1, varid_VertMode2, varid_VertMode3, varid_VertMode4, varid_VertMode5, varid_VertMode6
INTEGER                      :: varid_VertEV1, varid_VertEV2, varid_VertEV3, varid_VertEV4, varid_VertEV5, varid_VertEV6
INTEGER                      :: varid_HorizEV1, varid_HorizEV2, varid_HorizEV3, varid_HorizEV4, varid_HorizEV5, varid_HorizEV6
INTEGER                      :: varid_Regression, varid_Cov_rbalrbal, varid_Cov_rtotrbal

INTEGER                      :: ncid, ierr, i, x
INTEGER                      :: ierr1, ierr2, ierr3, ierr4, ierr5, ierr6, ierr7, ierr8, ierr9, ierr10, ierr11
INTEGER                      :: startA(1:1), countA(1:1)
INTEGER                      :: startB(1:2), countB(1:2)
INTEGER                      :: startC(1:3), countC(1:3)
REAL(ZREAL8)                 :: value



!*****************************************************************************************
!PRINT*, 'Read Covs'
!*****************************************************************************************


  ierr = NF_OPEN(filename, NF_NOWRITE, ncid)
  IF ( ierr .NE. 0 ) THEN
    PRINT*, ' *** Error opening file ***'
    PRINT*, ierr, NF_STRERROR(ierr)
    PRINT*, 'FILE :: ', filename
    STOP
  ENDIF


  !-----------------------------------------------------
  ! Input the values of the main variables as requested
  ! ----------------------------------------------------

  IF (input_params) THEN

    !--- Parameters and options ---
    ierr1   = NF_INQ_VARID(ncid, 'A', varid_A)
    ierr2   = NF_INQ_VARID(ncid, 'B', varid_B)
    ierr3   = NF_INQ_VARID(ncid, 'C', varid_C)
    ierr4   = NF_INQ_VARID(ncid, 'f', varid_f)
    ierr5   = NF_INQ_VARID(ncid, 'cvt_order', varid_CVT_order)
    ierr6   = NF_INQ_VARID(ncid, 'cvt_param_opt_gb', varid_CVT_param_opt_gb)
    ierr7   = NF_INQ_VARID(ncid, 'cvt_param_opt_hb', varid_CVT_param_opt_hb)
    ierr8   = NF_INQ_VARID(ncid, 'cvt_param_opt_ab', varid_CVT_param_opt_ab)
    ierr9   = NF_INQ_VARID(ncid, 'cvt_param_opt_reg', varid_CVT_param_opt_reg)
    ierr10  = NF_INQ_VARID(ncid, 'cvt_vert_opt_sym', varid_CVT_vert_opt_sym)
    ierr11  = NF_INQ_VARID(ncid, 'cvt_stddev_opt', varid_CVT_stddev_opt)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11) /= 0) THEN
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
      STOP
    END IF


    ierr1     = NF_GET_VAR1_DOUBLE(ncid, varid_A, 1, A)
    ierr2     = NF_GET_VAR1_DOUBLE(ncid, varid_B, 1, B)
    ierr3     = NF_GET_VAR1_DOUBLE(ncid, varid_C, 1, C)
    ierr4     = NF_GET_VAR1_DOUBLE(ncid, varid_f, 1, f)
    ierr5     = NF_GET_VAR1_INT(ncid, varid_CVT_order, 1, CVT % CVT_order)
    ierr6     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_gb, 1, CVT % CVT_param_opt_gb)
    ierr7     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_hb, 1, CVT % CVT_param_opt_hb)
    ierr8     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_ab, 1, CVT % CVT_param_opt_ab)
    ierr9     = NF_GET_VAR1_INT(ncid, varid_CVT_param_opt_reg, 1, CVT % CVT_param_opt_reg)
    ierr10    = NF_GET_VAR1_INT(ncid, varid_CVT_vert_opt_sym, 1, CVT % CVT_vert_opt_sym)
    ierr11    = NF_GET_VAR1_INT(ncid, varid_CVT_stddev_opt, 1, CVT % CVT_stddev_opt)
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 + ierr7 + ierr8 + ierr9 + ierr10 + ierr11) /= 0) THEN
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
      STOP
    END IF

    PRINT *, 'Read in A                 = ', A
    PRINT *, 'Read in B                 = ', B
    PRINT *, 'Read in C                 = ', C
    PRINT *, 'Read in f                 = ', f
    PRINT *, 'Read in CVT_order         = ', CVT % CVT_order
    PRINT *, 'Read in CVT_param_opt_gb  = ', CVT % CVT_param_opt_gb
    PRINT *, 'Read in CVT_param_opt_hb  = ', CVT % CVT_param_opt_hb
    PRINT *, 'Read in CVT_param_opt_ab  = ', CVT % CVT_param_opt_ab
    PRINT *, 'Read in CVT_param_opt_reg = ', CVT % CVT_param_opt_reg
    PRINT *, 'Read in CVT_vert_opt_sym  = ', CVT % CVT_vert_opt_sym
    PRINT *, 'Read in CVT_stddev_opt    = ', CVT % CVT_stddev_opt


  END IF



  IF (input_stddevs) THEN

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
      ierr1   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma1, 1, value)
      CVT % sigma1(1:nlongs,1:nlevs) = value
      ierr2   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma2, 1, value)
      CVT % sigma2(1:nlongs,1:nlevs) = value
      ierr3   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma3, 1, value)
      CVT % sigma3(1:nlongs,1:nlevs) = value
      ierr4   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma4, 1, value)
      CVT % sigma4(1:nlongs,1:nlevs) = value
      ierr5   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma5, 1, value)
      CVT % sigma5(1:nlongs,1:nlevs) = value
      ierr6   = NF_GET_VAR1_DOUBLE(ncid, varid_sigma6, 1, value)
      CVT % sigma6(1:nlongs,1:nlevs) = value

    CASE (2)
      ! Stddev level dependent
      startA(1) = 1
      countA(1) = nlevs
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_sigma1, startA, countA, CVT % sigma1(1, 1:nlevs))        
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_sigma2, startA, countA, CVT % sigma2(1, 1:nlevs))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_sigma3, startA, countA, CVT % sigma3(1, 1:nlevs))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_sigma4, startA, countA, CVT % sigma4(1, 1:nlevs))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_sigma5, startA, countA, CVT % sigma5(1, 1:nlevs))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_sigma6, startA, countA, CVT % sigma6(1, 1:nlevs))
      DO x = 2, nlongs
        CVT % sigma1(x, 1:nlevs) = CVT % sigma1(1, 1:nlevs)
        CVT % sigma2(x, 1:nlevs) = CVT % sigma2(1, 1:nlevs)
        CVT % sigma3(x, 1:nlevs) = CVT % sigma3(1, 1:nlevs)
        CVT % sigma4(x, 1:nlevs) = CVT % sigma4(1, 1:nlevs)
        CVT % sigma5(x, 1:nlevs) = CVT % sigma5(1, 1:nlevs)
        CVT % sigma6(x, 1:nlevs) = CVT % sigma6(1, 1:nlevs)
      END DO
    CASE (3)
      ! Stddev longitude and level dependent
      startB(1) = 1
      countB(1) = nlongs
      startB(2) = 1
      countB(2) = nlevs
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_sigma1, startB, countB, CVT % sigma1(1:nlongs, 1:nlevs))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_sigma2, startB, countB, CVT % sigma2(1:nlongs, 1:nlevs))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_sigma3, startB, countB, CVT % sigma3(1:nlongs, 1:nlevs))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_sigma4, startB, countB, CVT % sigma4(1:nlongs, 1:nlevs))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_sigma5, startB, countB, CVT % sigma5(1:nlongs, 1:nlevs))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_sigma6, startB, countB, CVT % sigma6(1:nlongs, 1:nlevs))

    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6) /= 0) THEN
      PRINT*, '***Error reading main variables (batch 2)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

  END IF


  IF (input_vert) THEN

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
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode1, startB, countB, CVT % VertMode1(1:nlevs, 1:nlevs, 1))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode2, startB, countB, CVT % VertMode2(1:nlevs, 1:nlevs, 1))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode3, startB, countB, CVT % VertMode3(1:nlevs, 1:nlevs, 1))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode4, startB, countB, CVT % VertMode4(1:nlevs, 1:nlevs, 1))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode5, startB, countB, CVT % VertMode5(1:nlevs, 1:nlevs, 1))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode6, startB, countB, CVT % VertMode6(1:nlevs, 1:nlevs, 1))

    CASE (2)
      ! Reversed MetO order of transforms (vert modes fn of horiz wavenumber)
      startC(1) = 1
      countC(1) = nlevs
      startC(2) = 1
      countC(2) = nlevs
      startC(3) = 1
      countC(3) = nlongs/2+1
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode1, startC, countC, CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode2, startC, countC, CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode3, startC, countC, CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode4, startC, countC, CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode5, startC, countC, CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_VertMode6, startC, countC, CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1))
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error reading main variables (batch 3)***'
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
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV1, startA, countA, CVT % VertEV1(1:nlevs, 1))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV2, startA, countA, CVT % VertEV2(1:nlevs, 1))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV3, startA, countA, CVT % VertEV3(1:nlevs, 1))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV4, startA, countA, CVT % VertEV4(1:nlevs, 1))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV5, startA, countA, CVT % VertEV5(1:nlevs, 1))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV6, startA, countA, CVT % VertEV6(1:nlevs, 1))

    CASE (2)
      ! Reversed MetO order of transforms (eigenvalues fn of horiz wavenumber)
      startB(1) = 1
      countB(1) = nlevs
      startB(2) = 1
      countB(2) = nlongs/2+1
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV1, startB, countB, CVT % VertEV1(1:nlevs, 1:nlongs/2+1))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV2, startB, countB, CVT % VertEV2(1:nlevs, 1:nlongs/2+1))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV3, startB, countB, CVT % VertEV3(1:nlevs, 1:nlongs/2+1))
      ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV4, startB, countB, CVT % VertEV4(1:nlevs, 1:nlongs/2+1))
      ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV5, startB, countB, CVT % VertEV5(1:nlevs, 1:nlongs/2+1))
      ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_VertEV6, startB, countB, CVT % VertEV6(1:nlevs, 1:nlongs/2+1))
    END SELECT

    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error reading main variables (batch 4)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF

  END IF


  IF (input_horiz) THEN
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
    ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV1, startB, countB, CVT % HorizEV1(1:nlongs/2+1, 1:nlevs))
    ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV2, startB, countB, CVT % HorizEV2(1:nlongs/2+1, 1:nlevs))
    ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV3, startB, countB, CVT % HorizEV3(1:nlongs/2+1, 1:nlevs))
    ierr4     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV4, startB, countB, CVT % HorizEV4(1:nlongs/2+1, 1:nlevs))
    ierr5     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV5, startB, countB, CVT % HorizEV5(1:nlongs/2+1, 1:nlevs))
    ierr6     = NF_GET_VARA_DOUBLE(ncid, varid_HorizEV6, startB, countB, CVT % HorizEV6(1:nlongs/2+1, 1:nlevs))
    IF ((ierr1 + ierr2 + ierr3 + ierr4 + ierr5 + ierr6 ) /= 0) THEN
      PRINT*, '***Error reading main variables (batch 5)***'
      PRINT*, 'ierr1 ',  ierr1,  NF_STRERROR(ierr1)
      PRINT*, 'ierr2 ',  ierr2,  NF_STRERROR(ierr2)
      PRINT*, 'ierr3 ',  ierr3,  NF_STRERROR(ierr3)
      PRINT*, 'ierr4 ',  ierr4,  NF_STRERROR(ierr4)
      PRINT*, 'ierr5 ',  ierr5,  NF_STRERROR(ierr5)
      PRINT*, 'ierr6 ',  ierr6,  NF_STRERROR(ierr6)
      STOP
    END IF
  END IF


  IF (input_regression) THEN

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
      ierr1     = NF_GET_VARA_DOUBLE(ncid, varid_Regression, startB, countB, CVT % Regression(1:nlevs, 1:nlevs))
      ierr2     = NF_GET_VARA_DOUBLE(ncid, varid_Cov_rbalrbal, startB, countB, CVT % Cov_rbalrbal(1:nlevs, 1:nlevs))
      ierr3     = NF_GET_VARA_DOUBLE(ncid, varid_Cov_rtotrbal, startB, countB, CVT % Cov_rtotrbal(1:nlevs, 1:nlevs))
      IF ((ierr1 + ierr2 + ierr3 ) /= 0) THEN
        PRINT*, '***Error reading main variables (batch 6)***'
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


END SUBROUTINE Read_Covs
