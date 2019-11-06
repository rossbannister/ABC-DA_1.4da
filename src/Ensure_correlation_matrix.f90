SUBROUTINE Ensure_correlation_matrix (L, mat)

! Make a covariance matrix into a correlation matrix
! (except when variances are zero)


USE DefConsTypes, ONLY :     &
    ZREAL8

IMPLICIT NONE

! Subroutine arguments
INTEGER,        INTENT(IN)      :: L              ! Size of matrix
REAL(ZREAL8),   INTENT(INOUT)   :: mat(1:L, 1:L)  ! in:  covariance matrix
                                                  ! out: correlation matrix
REAL(ZREAL8)                    :: srdiag(1:L)
INTEGER                         :: z1, z2


PRINT *, '*** Inside Ensure_correlation_matrix'
! Extract squareroot of diagonal elements
DO z1 = 1, L
  IF (mat(z1,z1) > 0.0) THEN
    srdiag(z1) = SQRT(mat(z1,z1))
  ELSE
    srdiag(z1) = 1.0
  END IF
END DO
PRINT *, 'sqrt(diag elements) (should be close to 1): ', srdiag(1:L)

! Normalise
DO z1 = 1, L
  DO z2 = 1, L
    mat(z1,z2) = mat(z1,z2) / (srdiag(z1) * srdiag(z2))
  END DO
END DO

END SUBROUTINE Ensure_correlation_matrix
