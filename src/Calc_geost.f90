SUBROUTINE Calc_geost (state)

! Calculate geostrophic imbalance in real space
USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  nlongs, nlevs,          &
  f, C, recipdx

IMPLICIT NONE

! Declare parameters
!-------------------
TYPE(ABC_type), INTENT(INOUT)    :: state

! Declare local variables
! -----------------------
INTEGER                          :: x, z
REAL(ZREAL8)                     :: rms1, rms2, norm
REAL(ZREAL8)                     :: term1(1:nlongs,1:nlevs), term2(1:nlongs,1:nlevs)

! Functions
! ---------
REAL(ZREAL8)                     :: RMS

! Calculate each term in the geostrophic imbalance equation
DO z = 1, nlevs
  DO x = 1, nlongs
    term1(x,z) = C * (state % r(x,z) - state % r(x-1,z)) * recipdx
    term2(x,z) = -1. * f * (state % v(x,z) + state % v(x-1,z)) / 2.
  END DO
END DO

! Find the RMS of each
rms1 = RMS(term1(1:nlongs,1:nlevs))
rms2 = RMS(term2(1:nlongs,1:nlevs))
norm = 1. / (rms1 + rms2)

! Compute the normalised geostrophic imbalance diagnostic
state % geost_imbal(1:nlongs,1:nlevs) = (term1(1:nlongs,1:nlevs) + term2(1:nlongs,1:nlevs)) * norm


END SUBROUTINE Calc_geost
