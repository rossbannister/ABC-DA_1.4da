SUBROUTINE Anbalw (LSrho, u, w_b, dims)

! Code to compute the anelastically balanced component of w
! Perturbations of density are neglected

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  recip2dx
  

IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)    :: LSrho(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8),    INTENT(IN)    :: u(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8),    INTENT(INOUT) :: w_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

REAL(ZREAL8)                   :: INT_HF
INTEGER                        :: x, z
REAL(ZREAL8)                   :: ux, uxm1, rhox, rhoxm1, rhoxp1, deltaz, drhodz
REAL(ZREAL8)                   :: diff(1:nlevs-1)
REAL(ZREAL8)                   :: rhow(0:nlevs)
REAL(ZREAL8)                   :: rho_full_levs(0:nlongs+1, 0:nlevs)


! Calculate LS rho on full levels
DO x = 0, nlongs+1
  DO z = 0, nlevs
    rho_full_levs(x,z) = INT_HF(LSrho(x,z), LSrho(x,z+1), z, dims)
  END DO
END DO


DO x = 1, nlongs
  ! Calculate the d(rho u)/dx for a column at this longitude
  DO z = 1, nlevs-1
    ! u interpolated to full-level
    ux      = INT_HF(u(x,z), u(x,z+1), z, dims)
    uxm1    = INT_HF(u(x-1,z), u(x-1,z+1), z, dims)
    ! rho interpolated to full-level
    rhox    = rho_full_levs(x,z)
    rhoxm1  = rho_full_levs(x-1,z)
    rhoxp1  = rho_full_levs(x+1,z)
    ! d(rho u)/dx
    diff(z) = recip2dx * ( (rhoxp1 + rhox) * ux - (rhox + rhoxm1) * uxm1 )
  END DO

  ! Integrate the continuity equation from the ground upwards
  rhow(0) = 0.0
  DO z = 1, nlevs-1
    ! Level width
    deltaz  = dims % half_levs(z+1) - dims % half_levs(z)
    rhow(z) = rhow(z-1) - diff(z) * deltaz
  END DO
  rhow(nlevs) = 0.0

  ! Compute balanced w
    w_b(x,1:nlevs) = rhow(1:nlevs) / rho_full_levs(x,1:nlevs)
END DO

END SUBROUTINE Anbalw
