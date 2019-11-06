SUBROUTINE Calc_vert_cov1 (Nmems, state1, state2, Cov)

! Calculate a vertical cov matrix for input data in format 1

USE DefConsTypes, ONLY :   &
    nlongs, nlevs,         &
    ZREAL8

IMPLICIT NONE

! Declare parameters
!---------------------
INTEGER,      INTENT(IN)    :: Nmems
REAL(ZREAL8), INTENT(IN)    :: state1(1:nlongs,1:nlevs, 1:Nmems)
REAL(ZREAL8), INTENT(IN)    :: state2(1:nlongs,1:nlevs, 1:Nmems)
REAL(ZREAL8), INTENT(INOUT) :: Cov(1:nlevs, 1:nlevs)

! Declare variables
!---------------------
INTEGER                     :: x, z1, z2, mem
REAL(ZREAL8)                :: total, RecipNmems, RecipNmemsm1
REAL(ZREAL8)                :: mean1(1:nlevs)
REAL(ZREAL8)                :: mean2(1:nlevs)

RecipNmems   = 1.0 / REAL(Nmems*nlongs)
RecipNmemsm1 = 1.0 / REAL(Nmems*nlongs-1)

! Compute the means
DO z1 = 1, nlevs
  mean1(z1) = SUM(state1(1:nlongs,z1,1:Nmems)) * RecipNmems
  mean2(z1) = SUM(state2(1:nlongs,z1,1:Nmems)) * RecipNmems
END DO

! Compute the covariances
DO z1 = 1, nlevs
  DO z2 = 1, nlevs
    total = 0.0
    DO mem = 1, Nmems
      DO x = 1, nlongs
        total = total + (state1(x,z1,mem) - mean1(z1)) *    &
                        (state2(x,z2,mem) - mean2(z2))
      END DO
    END DO
    Cov(z1,z2) = total * RecipNmemsm1
  ENDDO
ENDDO


END SUBROUTINE Calc_vert_cov1
