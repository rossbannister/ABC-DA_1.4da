SUBROUTINE Anbalw_adj (LSrho, u, w_b, dims)

! Code to compute the adjoint of the operator that computes the anelastically balanced component of w
! Perturbations of density are neglected

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  recip2dx
  

IMPLICIT NONE

REAL(ZREAL8),    INTENT(IN)    :: LSrho(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8),    INTENT(INOUT) :: u(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8),    INTENT(IN)    :: w_b(1:nlongs,1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

REAL(ZREAL8)                   :: INT_HF
INTEGER                        :: x, z
REAL(ZREAL8)                   :: ux, uxm1, rhox, rhoxm1, rhoxp1, deltaz, drhodz
REAL(ZREAL8)                   :: diff(1:nlevs)
REAL(ZREAL8)                   :: rhow(0:nlevs)
REAL(ZREAL8)                   :: rho_full_levs(0:nlongs+1, 0:nlevs)


! Calculate LS rho on full levels
DO x = 0, nlongs+1
  DO z = 0, nlevs
    rho_full_levs(x,z) = INT_HF(LSrho(x,z), LSrho(x,z+1), z, dims)
  END DO
END DO


DO x = 1, nlongs

  ! Compute balanced w
  rhow(1:nlevs) = w_b(x,1:nlevs) / rho_full_levs(x,1:nlevs)

  ! Integrate the continuity equation from the ground upwards
  diff(nlevs) = 0.0
  DO z = nlevs-1, 1, -1
    ! Level width
    deltaz  = dims % half_levs(z+1) - dims % half_levs(z)
    diff(z) = diff(z+1) - rhow(z) * deltaz
  END DO

  ! Calculate the d(rho u)/dx for a column at this longitude
  DO z = 1, nlevs-1
    ! rho interpolated to full-level
    rhox    = rho_full_levs(x,z)
    rhoxm1  = rho_full_levs(x-1,z)
    rhoxp1  = rho_full_levs(x+1,z)
    ! d(rho u)/dx
    ux   = recip2dx * (rhoxp1 + rhox) * diff(z)
    uxm1 = -1.0 * recip2dx * (rhox + rhoxm1) * diff(z)

    ! u interpolated to full-level
    CALL INT_HF_adj(u(x,z), u(x,z+1), ux, z, dims)
    !!ux      = INT_HF(u(x,z), u(x,z+1), z, dims)
    CALL INT_HF_adj(u(x-1,z), u(x-1,z+1), uxm1, z, dims)
    !!uxm1    = INT_HF(u(x-1,z), u(x-1,z+1), z, dims)
  END DO


END DO

END SUBROUTINE Anbalw_adj
