! R. Petrie, 3.0: 30-7-13
! R. Bannister  : 13-11-17

!=================================================================================================
 REAL(ZREAL8) FUNCTION STDEV (field)
!********************************************
!* Function to calculate standard deviation *
!* of a 2d field that has dimension         *
!* 1:nlongs, 1:nlevs                        *
!********************************************

USE DefConsTypes, ONLY  : &
  ZREAL8,                 &
  nlongs,                 &
  nlevs

IMPLICIT NONE

!Declare parameters
REAL(ZREAL8), INTENT(IN)          :: field(1:nlongs, 1:nlevs)
REAL(ZREAL8)                      :: ave, std, recip

recip = 1. / REAL(nlongs * nlevs)
ave   = recip * SUM(field(1:nlongs,1:nlevs))
std   = recip * SUM((field(1:nlongs,1:nlevs) - ave) * (field(1:nlongs,1:nlevs) - ave))

STDEV = sqrt(std)

END FUNCTION STDEV



!=================================================================================================
REAL(ZREAL8) FUNCTION GAUSS (std)
!********************************************
!* Function to calculate a gaussian         *
!* distributed random variable using the    *
!* Box-Mueller algorithm where the          *
!* distribution has a variance of std^2     *
!********************************************

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  pi

IMPLICIT NONE

!Declare parameters
REAL(ZREAL8), INTENT(IN)         :: std
REAL                             :: RAND
REAL(ZREAL8)                     :: u1, u2

u1 = RAND(0)
u2 = RAND(0)

gauss = sqrt(-2*log(u1))*cos(2*pi*u2)*std;

END FUNCTION GAUSS



!=================================================================================================
REAL(ZREAL8) FUNCTION RMS (state)

!********************************************
!* Function to calculate rms of a 2d-state  *
!* of dims (1:nlongs, 1:nlevs)              *
!********************************************

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  nlongs,                &
  nlevs


IMPLICIT NONE

! Declare parameters
!------------------
REAL(ZREAL8), INTENT(IN)    :: state(1:nlongs, 1:nlevs)

! Declare variables
!-------------------
REAL(ZREAL8)                :: eps, recip

eps   = 10E-19
recip = 1. / REAL(nlongs*nlevs)

! Calculate rmse
!-----------------
RMS = SQRT(recip * SUM(state(1:nlongs,1:nlevs) * state(1:nlongs,1:nlevs)))

IF (RMS .LT. eps) PRINT*, 'RMS SMALL'

END FUNCTION RMS




!=================================================================================================
REAL(ZREAL8) FUNCTION INT_FH (varl, varu, z, Dims)
! Function interpolates variable from full to half levs
! at a given grid point
! varl at (i,k-1)
! varu at (i,k)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  dims_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(IN) :: varl
REAL(ZREAL8),          INTENT(IN) :: varu
INTEGER,               INTENT(IN) :: z
TYPE(dims_type),       INTENT(IN) :: dims

! Calculate variable interpolated to full levs
INT_FH = Dims % a1(z) * varu + Dims % b1(z) * varl

END FUNCTION INT_FH



!!================================================================================================
REAL(ZREAL8) FUNCTION INT_HF (varl, varu, z, Dims)
! Function interpolates variable from half to full levs
! at a given grid point
! varl at (i,k)
! varu at (i,k+1)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  dims_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(IN) :: varl
REAL(ZREAL8),          INTENT(IN) :: varu
INTEGER,               INTENT(IN) :: z
TYPE(dims_type),       INTENT(IN) :: Dims

!Calculate variable interpolated to half levs
INT_HF = Dims % a2(z) * varu + Dims % b2(z) * varl

END FUNCTION INT_HF


!=================================================================================================
SUBROUTINE INT_FH_adj (varl, varu, interpolated, z, Dims)
! Function interpolates variable from full to half levs
! at a given grid point
! varl at (i,k-1)
! varu at (i,k)
! (Converted to a subroutine as there are two outputs from the adjoint routine)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  dims_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(INOUT) :: varl
REAL(ZREAL8),          INTENT(INOUT) :: varu
REAL(ZREAL8),          INTENT(IN)    :: interpolated
INTEGER,               INTENT(IN)    :: z
TYPE(dims_type),       INTENT(IN)    :: dims

! Calculate variable interpolated to full levs
varu = varu + Dims % a1(z) * interpolated
varl = varl + Dims % b1(z) * interpolated

END SUBROUTINE INT_FH_adj



!!================================================================================================
SUBROUTINE INT_HF_adj (varl, varu, interpolated, z, Dims)
! Function interpolates variable from half to full levs
! at a given grid point
! varl at (i,k)
! varu at (i,k+1)
! (Converted to a subroutine as there are two outputs from the adjoint routine)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  dims_type

IMPLICIT NONE

!Declare parameters
!------------------
REAL(ZREAL8),          INTENT(INOUT) :: varl
REAL(ZREAL8),          INTENT(INOUT) :: varu
REAL(ZREAL8),          INTENT(IN)    :: interpolated
INTEGER,               INTENT(IN)    :: z
TYPE(dims_type),       INTENT(IN)    :: Dims

!Calculate variable interpolated to half levs
varu = varu + Dims % a2(z) * interpolated
varl = varl + Dims % b2(z) * interpolated

END SUBROUTINE INT_HF_adj



!=================================================================================================
REAL(ZREAL8) FUNCTION InnerProdModelSpace (x1, x2, do_u, do_v, do_w, do_r, do_b, do_tracer, ignore_halos)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  ABC_type,              &
  nlongs,                &
  nlevs

IMPLICIT NONE

TYPE(ABC_type),     INTENT(IN) :: x1
TYPE(ABC_type),     INTENT(IN) :: x2
LOGICAL, OPTIONAL,  INTENT(IN) :: do_u
LOGICAL, OPTIONAL,  INTENT(IN) :: do_v
LOGICAL, OPTIONAL,  INTENT(IN) :: do_w
LOGICAL, OPTIONAL,  INTENT(IN) :: do_r
LOGICAL, OPTIONAL,  INTENT(IN) :: do_b
LOGICAL, OPTIONAL,  INTENT(IN) :: do_tracer
LOGICAL, OPTIONAL,  INTENT(IN) :: ignore_halos

REAL(ZREAL8)                   :: total_u, total_v, total_w, total_r, total_b, total_tracer
INTEGER                        :: xlow, xhigh, zlow, zhigh
LOGICAL                        :: flag_u, flag_v, flag_w, flag_r, flag_b, flag_tracer
LOGICAL                        :: some_flags, ignore_hs


IF (PRESENT(do_u)) THEN
  flag_u = do_u
ELSE
  flag_u = .FALSE.
END IF
IF (PRESENT(do_v)) THEN
  flag_v = do_v
ELSE
  flag_v = .FALSE.
END IF
IF (PRESENT(do_w)) THEN
  flag_w = do_w
ELSE
  flag_w = .FALSE.
END IF
IF (PRESENT(do_r)) THEN
  flag_r = do_r
ELSE
  flag_r = .FALSE.
END IF
IF (PRESENT(do_b)) THEN
  flag_b = do_b
ELSE
  flag_b = .FALSE.
END IF
IF (PRESENT(do_tracer)) THEN
  flag_tracer = do_tracer
ELSE
  flag_tracer = .FALSE.
END IF


some_flags = PRESENT(do_u) .OR. PRESENT(do_v) .OR. PRESENT(do_w) .OR. PRESENT(do_r) .OR. &
             PRESENT(do_b) .OR. PRESENT(do_tracer)

IF (.NOT.some_flags) THEN
  ! No flags set.  This is shorthand for all flags
  flag_u      = .TRUE.
  flag_v      = .TRUE.
  flag_w      = .TRUE.
  flag_r      = .TRUE.
  flag_b      = .TRUE.
  flag_tracer = .TRUE.
END IF

IF (PRESENT(ignore_halos)) THEN
  ignore_hs = ignore_halos
ELSE
  ignore_hs = .FALSE.
END IF

IF (ignore_hs) THEN
  xlow  = 1
  xhigh = nlongs
  zlow  = 1
  zhigh = nlevs
ELSE
  xlow  = 0
  xhigh = nlongs+1
  zlow  = 0
  zhigh = nlevs+1
END IF

IF (flag_u) THEN
  total_u      = SUM(x1 % u(xlow:xhigh,zlow:zhigh) * x2 % u(xlow:xhigh,zlow:zhigh))
ELSE
  total_u      = 0.0
END IF

IF (flag_v) THEN
  total_v      = SUM(x1 % v(xlow:xhigh,zlow:zhigh) * x2 % v(xlow:xhigh,zlow:zhigh))
ELSE
  total_v      = 0.0
END IF

IF (flag_w) THEN
  total_w      = SUM(x1 % w(xlow:xhigh,zlow:zhigh) * x2 % w(xlow:xhigh,zlow:zhigh))
ELSE
  total_w      = 0.0
END IF

IF (flag_r) THEN
  total_r      = SUM(x1 % r(xlow:xhigh,zlow:zhigh) * x2 % r(xlow:xhigh,zlow:zhigh))
ELSE
  total_r      = 0.0
END IF

IF (flag_b) THEN
  total_b      = SUM(x1 % b(xlow:xhigh,zlow:zhigh) * x2 % b(xlow:xhigh,zlow:zhigh))
ELSE
  total_b      = 0.0
END IF

IF (flag_tracer) THEN
  total_tracer = SUM(x1 % tracer(xlow:xhigh,zlow:zhigh) * x2 % tracer(xlow:xhigh,zlow:zhigh))
ELSE
  total_tracer = 0.0
END IF

InnerProdModelSpace = total_u + total_v + total_w + total_r + total_b + total_tracer

END FUNCTION InnerProdModelSpace



!=================================================================================================
COMPLEX(ZREAL8) FUNCTION InnerProdControlSpace (x1, x2, ComplexSpace, &
                                                do_1, do_2, do_3, do_4, do_5, do_6, ignore_halos)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  CV_type,               &
  nlongs,                &
  nlevs

IMPLICIT NONE

TYPE(CV_type),      INTENT(IN) :: x1
TYPE(CV_type),      INTENT(IN) :: x2
LOGICAL, OPTIONAL,  INTENT(IN) :: ComplexSpace
LOGICAL, OPTIONAL,  INTENT(IN) :: do_1
LOGICAL, OPTIONAL,  INTENT(IN) :: do_2
LOGICAL, OPTIONAL,  INTENT(IN) :: do_3
LOGICAL, OPTIONAL,  INTENT(IN) :: do_4
LOGICAL, OPTIONAL,  INTENT(IN) :: do_5
LOGICAL, OPTIONAL,  INTENT(IN) :: do_6
LOGICAL, OPTIONAL,  INTENT(IN) :: ignore_halos


COMPLEX(ZREAL8)                :: total_1, total_2, total_3, total_4, total_5, total_6
INTEGER                        :: xlow, xhigh, zlow, zhigh
LOGICAL                        :: flag_1, flag_2, flag_3, flag_4, flag_5, flag_6
LOGICAL                        :: some_flags, flag_Complex, ignore_hs

COMPLEX(ZREAL8)                :: InnerProdComplexSpace


IF (PRESENT(ComplexSpace)) THEN
  flag_Complex = ComplexSpace
ELSE
  flag_Complex = .FALSE.
END IF

IF (PRESENT(do_1)) THEN
  flag_1 = do_1
ELSE
  flag_1 = .FALSE.
END IF
IF (PRESENT(do_2)) THEN
  flag_2 = do_2
ELSE
  flag_2 = .FALSE.
END IF
IF (PRESENT(do_3)) THEN
  flag_3 = do_3
ELSE
  flag_3 = .FALSE.
END IF
IF (PRESENT(do_4)) THEN
  flag_4 = do_4
ELSE
  flag_4 = .FALSE.
END IF
IF (PRESENT(do_5)) THEN
  flag_5 = do_5
ELSE
  flag_5 = .FALSE.
END IF
IF (PRESENT(do_6)) THEN
  flag_6 = do_6
ELSE
  flag_6 = .FALSE.
END IF


some_flags = PRESENT(do_1) .OR. PRESENT(do_2) .OR. PRESENT(do_3) .OR. PRESENT(do_4) .OR. &
             PRESENT(do_5) .OR. PRESENT(do_6)

IF (.NOT.some_flags) THEN
  ! No flags set.  This is shorthand for all flags
  flag_1 = .TRUE.
  flag_2 = .TRUE.
  flag_3 = .TRUE.
  flag_4 = .TRUE.
  flag_5 = .TRUE.
  flag_6 = .TRUE.
END IF

IF (PRESENT(ignore_halos)) THEN
  ignore_hs = ignore_halos
ELSE
  ignore_hs = .FALSE.
END IF

IF (ignore_hs) THEN
  xlow  = 1
  xhigh = nlongs
  zlow  = 1
  zhigh = nlevs
ELSE
  xlow  = 0
  xhigh = nlongs+1
  zlow  = 0
  zhigh = nlevs+1
END IF


IF (flag_1) THEN
  IF (flag_Complex) THEN
    total_1 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v1(xlow:xhigh,zlow:zhigh), x2 % v1(xlow:xhigh,zlow:zhigh))
  ELSE
    total_1 = COMPLEX(SUM(x1 % v1(xlow:xhigh,zlow:zhigh) * x2 % v1(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_1 = COMPLEX(0.0, 0.0)
END IF

IF (flag_2) THEN
  IF (flag_Complex) THEN
    total_2 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v2(xlow:xhigh,zlow:zhigh), x2 % v2(xlow:xhigh,zlow:zhigh))
  ELSE
    total_2 = COMPLEX(SUM(x1 % v2(xlow:xhigh,zlow:zhigh) * x2 % v2(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_2 = COMPLEX(0.0, 0.0)
END IF

IF (flag_3) THEN
  IF (flag_Complex) THEN
    total_3 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v3(xlow:xhigh,zlow:zhigh), x2 % v3(xlow:xhigh,zlow:zhigh))
  ELSE
    total_3 = COMPLEX(SUM(x1 % v3(xlow:xhigh,zlow:zhigh) * x2 % v3(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_3 = COMPLEX(0.0, 0.0)
END IF

IF (flag_4) THEN
  IF (flag_Complex) THEN
    total_4 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v4(xlow:xhigh,zlow:zhigh), x2 % v4(xlow:xhigh,zlow:zhigh))
  ELSE
    total_4 = COMPLEX(SUM(x1 % v4(xlow:xhigh,zlow:zhigh) * x2 % v4(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_4 = COMPLEX(0.0, 0.0)
END IF

IF (flag_5) THEN
  IF (flag_Complex) THEN
    total_5 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v5(xlow:xhigh,zlow:zhigh), x2 % v5(xlow:xhigh,zlow:zhigh))
  ELSE
    total_5 = COMPLEX(SUM(x1 % v5(xlow:xhigh,zlow:zhigh) * x2 % v5(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_5 = COMPLEX(0.0, 0.0)
END IF

IF (flag_6) THEN
  IF (flag_Complex) THEN
    total_6 = InnerProdComplexSpace(xlow, xhigh, zlow, zhigh, x1 % v6(xlow:xhigh,zlow:zhigh), x2 % v6(xlow:xhigh,zlow:zhigh))
  ELSE
    total_6 = COMPLEX(SUM(x1 % v6(xlow:xhigh,zlow:zhigh) * x2 % v6(xlow:xhigh,zlow:zhigh)), 0.0)
  END IF
ELSE
  total_6 = COMPLEX(0.0, 0.0)
END IF

InnerProdControlSpace = total_1 + total_2 + total_3 + total_4 + total_5 + total_6

END FUNCTION InnerProdControlSpace



!=================================================================================================
COMPLEX(ZREAL8) FUNCTION InnerProdComplexSpace (xlow, xhigh, zlow, zhigh, x1, x2)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  nlongs,                &
  nlevs

IMPLICIT NONE

INTEGER,      INTENT(IN) :: xlow
INTEGER,      INTENT(IN) :: xhigh
INTEGER,      INTENT(IN) :: zlow
INTEGER,      INTENT(IN) :: zhigh
REAL(ZREAL8), INTENT(IN) :: x1(xlow:xhigh,zlow:zhigh)
REAL(ZREAL8), INTENT(IN) :: x2(xlow:xhigh,zlow:zhigh)

INTEGER                  :: k, real_index, imag_index
COMPLEX(ZREAL8)          :: total

! Deal with the largest scale (this is real)
real_index = 1
total = COMPLEX(SUM(x1(real_index,zlow:zhigh) * x2(real_index,zlow:zhigh)), 0.0)
!total = COMPLEX(x1(real_index,1) * x2(real_index,1), 0.0)
! Deal with the bulk of the scales
DO k = 1, nlongs/2-1
  real_index = 2*k
  imag_index = 2*k+1

! Contribution from positive wavenumbers
!  total      = total + COMPLEX(SUM(x1(real_index,zlow:zhigh)*x2(real_index,zlow:zhigh)) +   &
!                               SUM(x1(imag_index,zlow:zhigh)*x2(imag_index,zlow:zhigh)),    &
!                               SUM(x1(real_index,zlow:zhigh)*x2(imag_index,zlow:zhigh)) -   &
!                               SUM(x1(imag_index,zlow:zhigh)*x2(real_index,zlow:zhigh))) / 4.
! Contribution from negative wavenumbers
!  total      = total + COMPLEX(SUM(x1(real_index,zlow:zhigh)*x2(real_index,zlow:zhigh)) +   &
!                               SUM(x1(imag_index,zlow:zhigh)*x2(imag_index,zlow:zhigh)),    &
!                               SUM(x1(imag_index,zlow:zhigh)*x2(real_index,zlow:zhigh)) -   &
!                               SUM(x1(real_index,zlow:zhigh)*x2(imag_index,zlow:zhigh))) / 4.

  ! Can replace the above commented-out code with the following
  ! Contribution from positive and negative wavenumbers in one (imag parts of above cancel)
  total      = total + COMPLEX(SUM(x1(real_index,zlow:zhigh)*x2(real_index,zlow:zhigh)) +   &
                               SUM(x1(imag_index,zlow:zhigh)*x2(imag_index,zlow:zhigh)),    &
                               0.0) / 2.



END DO
! Deal with the smallest scale
real_index = nlongs
total = total + COMPLEX(SUM(x1(real_index,zlow:zhigh) * x2(real_index,zlow:zhigh)), 0.0)

InnerProdComplexSpace = total

END FUNCTION InnerProdComplexSpace



!=================================================================================================
REAL(ZREAL8) FUNCTION InnerProduct_obsself ( Obs, weighted, part )
USE DefConsTypes, ONLY : &
  ZREAL8,     &
  Obs_type

! Computes the inner product of Obs (with itself)

IMPLICIT NONE

TYPE(Obs_type), POINTER, INTENT(IN) :: Obs
LOGICAL,                 INTENT(IN) :: weighted   ! .TRUE. if inv variance weighted
INTEGER,                 INTENT(IN) :: part       ! 1 means use deltay_m
                                                  ! 2 means use hxmy

! Note for adjoint tests, use weight=.FALSE., part=1
!      for Jo calc,       use weight=.TRUE., part=2

REAL(ZREAL8)                        :: total
TYPE(Obs_type), POINTER             :: thisob

total  = 0.0
thisob => Obs
DO
  IF (ASSOCIATED(thisob)) THEN
    IF (thisob % ob_ok) THEN
      ! This indicates that the tstep_lower, xbox_lower, and zbox_lower are already computed
      IF (weighted) THEN
        IF (part == 1) THEN
          total = total + thisob % deltay_m * thisob % deltay_m / thisob % variance
        ELSE
          total = total + thisob % hxmy * thisob % hxmy / thisob % variance
        END IF
      ELSE
        IF (part == 1) THEN
          total = total + thisob % deltay_m * thisob % deltay_m
        ELSE
          total = total + thisob % hxmy * thisob % hxmy
        END IF
      END IF
    END IF
    thisob => thisob % next
  ELSE
    EXIT
  END IF
END DO

InnerProduct_obsself = total

END FUNCTION InnerProduct_obsself



!=================================================================================================
INTEGER FUNCTION FindLowerIndex (value, arraysize, array)

USE DefConsTypes, ONLY : &
  ZREAL8

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN) :: value
INTEGER,      INTENT(IN) :: arraysize
REAL(ZREAL8), INTENT(IN) :: array(0:arraysize)

INTEGER                  :: testpoint


testpoint = -1
DO
  testpoint = testpoint + 1
  IF (testpoint .GE. arraysize) THEN
    testpoint = -1
    EXIT
  ELSE
    IF (array(testpoint) .GT. value) THEN
      ! This is the first point in the array that is greater than the value
      testpoint = testpoint - 1
      EXIT
    END IF
  END IF
END DO

FindLowerIndex = testpoint

END FUNCTION FindLowerIndex



!=================================================================================================
REAL(ZREAL8) FUNCTION Interpolate3D (tstates, States, xax, zax, tax, &
                                     xbox_lower, zbox_lower,         &
                                     xval, zval, tval,               &
                                     ob_of_what)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  ABC_type

IMPLICIT NONE

INTEGER,        INTENT(IN) :: tstates                ! The num of time states for interp (1 or 2)
TYPE(ABC_type), INTENT(IN) :: States(1:tstates)      ! The main field to be interpolated
REAL(ZREAL8),   INTENT(IN) :: xax(0:1)               ! The x-axis values
REAL(ZREAL8),   INTENT(IN) :: zax(0:1)               ! The z-axis values
INTEGER,        INTENT(IN) :: tax(1:tstates)         ! The t-axis values
INTEGER,        INTENT(IN) :: xbox_lower             ! Lower x-box index
INTEGER,        INTENT(IN) :: zbox_lower             ! Lower z-box index
REAL(ZREAL8),   INTENT(IN) :: xval                   ! The x point to interpolate to
REAL(ZREAL8),   INTENT(IN) :: zval                   ! The z point to interpolate to
INTEGER,        INTENT(IN) :: tval                   ! The t point to interpolate to
INTEGER,        INTENT(IN) :: ob_of_what             ! The field kind

REAL(ZREAL8)               :: Field2d(0:1,0:1)
REAL(ZREAL8)               :: Field1d(0:2)
INTEGER                    :: t, x, z
REAL(ZREAL8)               :: tax_real(1:tstates), tval_real

REAL(ZREAL8)               :: Interpolate1D


tax_real(1:tstates) = REAL(tax(1:tstates))
tval_real           = REAL(tval)

! Interpolate in time
SELECT CASE (ob_of_what)
CASE (1) ! u
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % u(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
      END IF
    END DO
  END DO

CASE (2) ! v
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % v(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
      END IF
    END DO
  END DO

CASE (3) ! w
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % w(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
      END IF
    END DO
  END DO

CASE (4) ! r
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % r(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
      END IF
    END DO
  END DO

CASE (5) ! b
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % b(xbox_lower+x, zbox_lower+z)
        !print *, '====='
        !PRINT *, xbox_lower+x, zbox_lower+z, States(t) % b(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        !PRINT *, x, z, ' | ', Field1d(1:2), ' | ', tax_real(1:2), ' | ', tval_real, ' | ', Field2d(x,z)
      END IF
    END DO
  END DO

CASE (6) ! tracer
  DO x = 0, 1
    DO z = 0, 1
      DO t = 1, tstates
        Field1d(t)   = States(t) % tracer(xbox_lower+x, zbox_lower+z)
      END DO
      IF (tstates == 1) THEN
        Field2d(x,z) = Field1d(1)
      ELSE
        Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
      END IF
    END DO
  END DO

END SELECT

! Interpolate in z
DO x = 0, 1
  Field1d(x) = Interpolate1D (Field2d(x,0:1), zax(0:1), zval)
END DO

! Interpolate in x
Interpolate3D = Interpolate1D (Field1d(0:1), xax(0:1), xval)

END FUNCTION Interpolate3D


!=================================================================================================
SUBROUTINE Interpolate3D_adj (ModelOb,                        &
                              tstates, States, xax, zax, tax, &
                              xbox_lower, zbox_lower,         &
                              xval, zval, tval,               &
                              ob_of_what)

! Adjoint of the function Interpolate3D
! ModelOb was the output of the forward function (input to adjoint)

USE DefConsTypes, ONLY : &
  ZREAL8,                &
  ABC_type,              &
  zero

IMPLICIT NONE

REAL(ZREAL8),   INTENT(IN)    :: ModelOb                ! Model Observation
INTEGER,        INTENT(IN)    :: tstates                ! The num of time states for interp (1 or 2)
TYPE(ABC_type), INTENT(INOUT) :: States(1:tstates)      ! The main field to be interpolated
REAL(ZREAL8),   INTENT(IN)    :: xax(0:1)               ! The x-axis values
REAL(ZREAL8),   INTENT(IN)    :: zax(0:1)               ! The z-axis values
INTEGER,        INTENT(IN)    :: tax(1:tstates)         ! The t-axis values
INTEGER,        INTENT(IN)    :: xbox_lower             ! Lower x-box index
INTEGER,        INTENT(IN)    :: zbox_lower             ! Lower z-box index
REAL(ZREAL8),   INTENT(IN)    :: xval                   ! The x point to interpolate to
REAL(ZREAL8),   INTENT(IN)    :: zval                   ! The z point to interpolate to
INTEGER,        INTENT(IN)    :: tval                   ! The t point to interpolate to
INTEGER,        INTENT(IN)    :: ob_of_what             ! The field kind

REAL(ZREAL8)                  :: Field2d(0:1,0:1)
REAL(ZREAL8)                  :: Field1d(0:2)
REAL(ZREAL8)                  :: Interpolate1D
INTEGER                       :: t, x, z
REAL(ZREAL8)                  :: tax_real(1:tstates), tval_real


tax_real(1:tstates) = REAL(tax(1:tstates))
tval_real           = REAL(tval)


! Forward parts are commented out

! Interpolate in x
!ModelOb = Interpolate1D (Field1d(0:1), xax(0:1), xval)

Field1d(0:2) = zero
CALL Interpolate1D_adj (ModelOb, Field1d(0:1), xax(0:1), xval)

! Interpolate in z
DO x = 0, 1
  !Field1d(x) = Interpolate1D (Field2d(x,0:1), zax(0:1), zval)

  Field2d(x,0:1) = zero
  CALL Interpolate1D_adj (Field1d(x), Field2d(x,0:1), zax(0:1), zval)
END DO


! Interpolate in time
SELECT CASE (ob_of_what)
CASE (1) ! u
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
         CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % u(xbox_lower+x, zbox_lower+z)
        States(t) % u(xbox_lower+x, zbox_lower+z) = States(t) % u(xbox_lower+x, zbox_lower+z) + &
                                                    Field1d(t)
      END DO
    END DO
  END DO

CASE (2) ! v
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % v(xbox_lower+x, zbox_lower+z)
        States(t) % v(xbox_lower+x, zbox_lower+z) = States(t) % v(xbox_lower+x, zbox_lower+z) + &
                                                    Field1d(t)
      END DO
    END DO
  END DO

CASE (3) ! w
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % w(xbox_lower+x, zbox_lower+z)
        States(t) % w(xbox_lower+x, zbox_lower+z) = States(t) % w(xbox_lower+x, zbox_lower+z) + &
                                                    Field1d(t)
      END DO
    END DO
  END DO

CASE (4) ! r
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % r(xbox_lower+x, zbox_lower+z)
        States(t) % r(xbox_lower+x, zbox_lower+z) = States(t) % r(xbox_lower+x, zbox_lower+z) + &
                                                    Field1d(t)
      END DO
    END DO
  END DO

CASE (5) ! b
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % b(xbox_lower+x, zbox_lower+z)
        States(t) % b(xbox_lower+x, zbox_lower+z) = States(t) % b(xbox_lower+x, zbox_lower+z) + &
                                                    Field1d(t)
      END DO
    END DO
  END DO

CASE (6) ! tracer
  DO x = 0, 1
    DO z = 0, 1
      Field1d(1:2) = zero
      IF (tstates == 1) THEN
        !Field2d(x,z) = Field1d(1)
        Field1d(1) = Field2d(x,z)
      ELSE
        !Field2d(x,z) = Interpolate1D (Field1d(1:2), tax_real(1:2), tval_real)
        CALL Interpolate1D_adj (Field2d(x,z), Field1d(1:2), tax_real(1:2), tval_real)
      END IF
      DO t = 1, tstates
        !Field1d(t)   = States(t) % tracer(xbox_lower+x, zbox_lower+z)
        States(t) % tracer(xbox_lower+x, zbox_lower+z) = States(t) % tracer(xbox_lower+x, zbox_lower+z) + &
                                                         Field1d(t)
      END DO
    END DO
  END DO

END SELECT

END SUBROUTINE Interpolate3D_adj



!=================================================================================================
REAL(ZREAL8) FUNCTION Interpolate1D (Field, axis, value)

USE DefConsTypes, ONLY : &
  ZREAL8

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN) :: Field(1:2)   ! The main field to be interpolated
REAL(ZREAL8), INTENT(IN) :: axis(1:2)    ! The axis values
REAL(ZREAL8), INTENT(IN) :: value        ! The point to interpolate to

REAL(ZREAL8)             :: dx, dy, dxpoint

dx      = axis(2) - axis(1)
dy      = Field(2) - Field(1)
dxpoint = value - axis(1)

Interpolate1D = Field(1) + dxpoint * dy / dx

END FUNCTION Interpolate1D




!=================================================================================================
SUBROUTINE Interpolate1D_adj (ModelOb, Field, axis, value)

! Adjoint of the function Interpolate1D
! ModelOb was the output of the forward function (input to adjoint)


USE DefConsTypes, ONLY : &
  ZREAL8

IMPLICIT NONE

REAL(ZREAL8), INTENT(IN)    :: ModelOb      ! Model Observation
REAL(ZREAL8), INTENT(INOUT) :: Field(1:2)   ! The main field to be interpolated
REAL(ZREAL8), INTENT(IN)    :: axis(1:2)    ! The axis values
REAL(ZREAL8), INTENT(IN)    :: value        ! The point to interpolate to

REAL(ZREAL8)                :: dx, dy, dxpoint

dx      = axis(2) - axis(1)
dxpoint = value - axis(1)

! Forward code commented out
!dy      = Field(2) - Field(1)
!ModelOb = Field(1) + dxpoint * dy / dx

! Adjoint of the above
Field(1) = Field(1) + ModelOb
dy       = ModelOb * dxpoint / dx
Field(2) = Field(2) + dy
Field(1) = Field(1) - dy


END SUBROUTINE Interpolate1D_adj
