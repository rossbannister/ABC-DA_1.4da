SUBROUTINE Calc_hydro (state, dims)

! Calculate hydrostatic imbalance in real space

USE DefConsTypes, ONLY :   &
    ABC_type,              &
    dims_type,             &
    ZREAL8,                &
    nlongs, nlevs, C

IMPLICIT NONE

! Declare parameters
!---------------------
TYPE(ABC_type),  INTENT(INOUT) :: state
TYPE(dims_type), INTENT(IN)    :: dims

! Declare variables
!---------------------
INTEGER                        :: x, z
REAL(ZREAL8)                   :: rms1, rms2, norm
REAL(ZREAL8)                   :: term1(1:nlongs,1:nlevs), term2(1:nlongs,1:nlevs)

! Functions
! ---------
REAL(ZREAL8)                   :: RMS

! Calculate each term in the hydrostatic imbalance equation
DO x = 1, nlongs
  DO z = 1, nlevs
    term1(x,z) = C * (state % r(x,z+1) - state % r(x,z)) /         &
                     (dims % half_levs(z+1) - dims % half_levs(z))
    term2(x,z) = -1. * state % b(x,z)
  ENDDO
ENDDO

! Find the RMS of each
rms1 = RMS(term1(1:nlongs,1:nlevs))
rms2 = RMS(term2(1:nlongs,1:nlevs))
norm = 1. / (rms1 + rms2)

! Compute the normalised geostrophic imbalance diagnostic
state % hydro_imbal(1:nlongs,1:nlevs) = (term1(1:nlongs,1:nlevs) + term2(1:nlongs,1:nlevs)) * norm

END SUBROUTINE Calc_hydro
