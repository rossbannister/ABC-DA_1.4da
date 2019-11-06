SUBROUTINE Energy (state)

! To calculate components of energy

USE DefConsTypes, ONLY :  &
    ABC_type,             &
    ZREAL8,               &
    nlongs, nlevs,        &
    A, B, C

IMPLICIT NONE

TYPE(ABC_type),  INTENT(INOUT)   :: state
INTEGER                          :: x, z
REAL(ZREAL8)                     :: u, v, w, rho, bp, r
REAL(ZREAL8)                     :: E_k, E_b, E_e

! Initialise
E_k = 0.0
E_e = 0.0
E_b = 0.0

! Calculate energy
DO z = 1, nlevs
  DO x = 1, nlongs
    r    = state % r(x,z)
    u    = state % u(x,z)
    v    = state % v(x,z)
    w    = state % w(x,z)
    bp   = state % b(x,z)
    rho  = state % rho(x,z)

    E_k  = E_k + rho * (u*u + v*v + w*w) / 2.0
    E_e  = E_e + C * r*r / (2. * B)
    E_b  = E_b + rho * bp*bp / (2. * A*A)
  END DO
END DO

state % Kinetic_Energy = E_k
state % Buoyant_Energy = E_b
state % Elastic_Energy = E_e
state % Total_Energy   = E_k + E_e + E_b

END SUBROUTINE Energy
