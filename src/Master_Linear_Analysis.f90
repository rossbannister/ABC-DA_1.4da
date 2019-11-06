PROGRAM Master_Linear_Analysis

!*****************************************************
!*   Code to calculate the eigenvalues and           *
!*   speeds fo the waves in the model and to         *
!*   detect the sensitivity of these to model        *
!*   parameters                                      *
!*                                                   *
!*   R. Petrie, R. Bannister                         *
!*   Version 3                                       *
!*   09-11-7                                        *
!*   Now using time units with dimensions            *
!*   (not scaled by f) 08-12-16                      *
!*****************************************************

USE DefConsTypes, ONLY :  &
  datadirLinearAnal,      &
  ZREAL8,                 &
  nlongs,                 &
  nlevs,                  &
  Lx,                     &
  pi,                     &
  H,                      &
  A,                      &
  B,                      &
  C,                      &
  f

IMPLICIT NONE


! DECLARE VARIABLES
!------------------
REAL(ZREAL8)    :: inv_L, inv_H, sqrtBC, kfact, mfact    
REAL(ZREAL8)    :: mall(0:nlevs-1), kall(0:nlongs-1)     
REAL(ZREAL8)    :: L(1:5, 1:5)                           
REAL(ZREAL8)    :: lambda(1:5)                           
REAL(ZREAL8)    :: Work(1:14)
REAL(ZREAL8)    :: grav_frequency (0:nlongs-1, 0:nlevs-1)    
REAL(ZREAL8)    :: acou_frequency (0:nlongs-1, 0:nlevs-1)    
REAL(ZREAL8)    :: hori_grav_speed(0:nlongs-1, 0:nlevs-1)    
REAL(ZREAL8)    :: hori_acou_speed(0:nlongs-1, 0:nlevs-1)    
REAL(ZREAL8)    :: vert_grav_speed(0:nlongs-1, 0:nlevs-1)    
REAL(ZREAL8)    :: vert_acou_speed(0:nlongs-1, 0:nlevs-1)    
INTEGER         :: k, m, i, j
INTEGER         :: ierr

PRINT*, '*************************************************************************'
PRINT*, 'Running Master_Linear_Analysis'
PRINT*, '*************************************************************************'

! Read namelist
CALL SetOptions

PRINT*, 'Lx = ', Lx

inv_L   = 1.0 / Lx
inv_H   = 1.0 / H
sqrtBC  = SQRT(B * C)
kfact   = Lx / (2. * pi)
mfact   = H / (2. * pi)


! Allocate horizontal and vertical wavenumbers
!---------------------------------------------
DO k = 0, nlongs-1
  kall(k) = 2. * REAL(k) * pi
END DO

DO m = 0, nlevs-1
  mall(m) = 2. * REAL(m) * pi
END DO



! Elements of matrix are computed below
DO m = 0, nlevs-1

  DO k = 0, nlongs-1

    L(1,1) = 0.0
    L(1,2) = f
    L(1,3) = 0.0
    L(1,4) = -1. * kall(k) * sqrtBC * inv_L
    L(1,5) = 0.0

    L(2,1) = L(1,2)
    L(2,2) = 0.0
    L(2,3) = 0.0
    L(2,4) = 0.0
    L(2,5) = 0.0

    L(3,1) = 0.0
    L(3,2) = 0.0
    L(3,3) = 0.0
    L(3,4) = mall(m) * sqrtBC * inv_H
    L(3,5) = A

    L(4,1) = L(1,4)
    L(4,2) = 0.0
    L(4,3) = L(3,4)
    L(4,4) = 0.0
    L(4,5) = 0.0

    L(5,1) = 0.0
    L(5,2) = 0.0
    L(5,3) = L(3,5)
    L(5,4) = 0.0
    L(5,5) = 0.0

    IF ( ((m == 0) .AND. (k == 0)) .OR. ((m == 3) .AND. (k == 7)) ) THEN
      PRINT *, 'Example matrix'
      DO i = 1, 5
        PRINT '(5F8.3)', (L(i,j), j=1,5)
      END DO
    END IF

    CALL DSYEV('N',         & ! Eigenvalues only to be computed
               'U',         & ! Upper triangular of matrix specified
               5,           & ! Order of the matrix to be diagonalised
               L(1:5, 1:5), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
               5,           & ! Leading dimension of the matrix
               lambda(1:5), & ! Eigenvalue array
               Work(1:14),  & ! Work array
               14,          & ! Size of work array
               ierr)
    IF (ierr /= 0) THEN
      PRINT*, 'Error with eigen routine', m, k
      STOP
    END IF

    acou_frequency(k,m) = lambda(5)  ! eig routine lists frequencies in order - - 0 + +
    grav_frequency(k,m) = lambda(4)  ! eig routine lists frequencies in order - - 0 + +

  END DO
END DO


! Calculate wave speeds
!----------------------

DO m = 0, nlevs-1
  DO k = 1, nlongs-1
    hori_grav_speed(k,m) = (grav_frequency(k,m) - grav_frequency(k-1,m))  * kfact
    hori_acou_speed(k,m) = (acou_frequency(k,m) - acou_frequency(k-1,m))  * kfact
  ENDDO
  hori_grav_speed(0,m)   = 0.0
  hori_acou_speed(0,m)   = 0.0
ENDDO

DO k = 0, nlongs-1
  DO m = 1, nlevs-1
    vert_grav_speed(k,m) = (grav_frequency(k,m)  - grav_frequency(k,m-1))  * mfact
    vert_acou_speed(k,m) = (acou_frequency(k,m)  - acou_frequency(k,m-1))  * mfact
  ENDDO
  vert_grav_speed(k,0)   = 0.0
  vert_acou_speed(k,0)   = 0.0

ENDDO

! Dump data
!------------


OPEN (51, file = TRIM(datadirLinearAnal) // '/grav_frequency.dat')
OPEN (52, file = TRIM(datadirLinearAnal) // '/acou_frequency.dat')
OPEN (53, file = TRIM(datadirLinearAnal) // '/hori_grav_speed.dat')
OPEN (54, file = TRIM(datadirLinearAnal) // '/hori_acou_speed.dat')
OPEN (55, file = TRIM(datadirLinearAnal) // '/vert_grav_speed.dat')
OPEN (56, file = TRIM(datadirLinearAnal) // '/vert_acou_speed.dat')


DO k = 0, nlongs-1
  WRITE (51, *) grav_frequency(k, 0:nlevs-1)
  WRITE (52, *) acou_frequency(k, 0:nlevs-1)
ENDDO

DO m = 0, nlevs-1
  WRITE (53, *) hori_grav_speed(0:nlongs-1, m)
  WRITE (54, *) hori_acou_speed(0:nlongs-1, m)
ENDDO

DO k = 0, nlongs-1
  WRITE (55, *) vert_grav_speed(k, 0:nlevs-1)
  WRITE (56, *) vert_acou_speed(k, 0:nlevs-1)
ENDDO


CLOSE(51)
CLOSE(52)
CLOSE(53)
CLOSE(54)
CLOSE(55)
CLOSE(56)

END PROGRAM Master_Linear_Analysis
