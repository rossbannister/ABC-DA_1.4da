SUBROUTINE InverseSymMat (dim, mat)

! Compute the inverse of a symmetrix matrix (with some conditioning adjustment)
! Uses eigendecomposition method

USE DefConsTypes, ONLY :   &
    ZREAL8,                &
    ConditionFudge

IMPLICIT NONE

! Declare parameters
!---------------------
INTEGER,      INTENT(IN)    :: dim
REAL(ZREAL8), INTENT(INOUT) :: mat(1:dim,1:dim)

! Declare variables
!---------------------
INTEGER                     :: ierr, i, j
REAL(ZREAL8)                :: evecs(1:dim,1:dim)
REAL(ZREAL8)                :: evals(1:dim)
REAL(ZREAL8)                :: Work(1:3*dim)
REAL(ZREAL8)                :: totaleigs, add, inv
REAL(ZREAL8)                :: invlambdaFT(1:dim,1:dim)



! Initial state for evecs is the matrix itself
evecs(1:dim,1:dim) = mat(1:dim,1:dim)

! Compute the eigenvalues and eigenvectors of the input matrix
CALL DSYEV('V',                 & ! Eigenvalues and vectors to be computed
           'U',                 & ! Upper triangular of matrix specified
           dim,                 & ! Order of the matrix to be diagonalised
           evecs(1:dim, 1:dim), & ! IN - matrix to be diagonalised, OUT - the eigenvectors
           dim,                 & ! Leading dimension of the matrix
           evals(1:dim),        & ! Eigenvalue array
           Work(1:3*dim),       & ! Work array
           3*dim,               & ! Size of work array
           ierr)

IF (ierr /= 0) THEN
  PRINT *, 'Error computing eigenvalues in InverseSymMat.f90'
  STOP
END IF

!PRINT *, 'Here is the eigenvale spectrum of the balanced r cov matrix'
!PRINT *, evals(1:dim)

! Condition the eigenvalues
totaleigs    = SUM(evals(1:dim))
add          = totaleigs * ConditionFudge
evals(1:dim) = evals(1:dim) + add

!PRINT *, 'Here is the modified eigenvale spectrum of the balanced r cov matrix'
!PRINT *, evals(1:dim)

DO i = 1, dim
  inv = 1.0 / evals(i)
  DO j = 1, dim
    invlambdaFT(i,j) = inv * evecs(j,i)
  END DO
END DO

mat(1:dim,1:dim) = MATMUL(evecs(1:dim,1:dim), invlambdaFT(1:dim,1:dim))


END SUBROUTINE InverseSymMat
