SUBROUTINE Effective_buoyancy (state, dims)

!**************************************************
!* Subroutine to calculate the effective buoyancy *
!**************************************************

USE DefConsTypes, ONLY :   &
    ABC_type,              &
    nlongs, nlevs,         &
    dims_type,             &
    A

IMPLICIT NONE

TYPE(ABC_type),  INTENT(INOUT)   :: state
TYPE(dims_type), INTENT(IN)      :: dims

! Declare local variables
!-------------------
INTEGER                          :: i, k

! Routine written assuming constant dz
! modification required if using Charney-Philips

! NOT NORMALIZED

DO k = 1, nlevs
  DO i = 1, nlongs

   state % b_ef(i,k) = (state % b(i,k+1) - state % b(i,k-1)) /           &
                       (dims % full_levs(k+1) - dims % full_levs(k-1)) + &
                       A*A

  END DO
END DO

END SUBROUTINE Effective_buoyancy
