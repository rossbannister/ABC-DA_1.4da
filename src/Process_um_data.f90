SUBROUTINE Process_um_data (um_data, state, dims)

!************************************
!* Subroutine to process um         *
!* data for toy model               *
!*                                  *
!* Version history                  *
!* R. Petrie, 1.0: 6-6-2011         *
!* R. Petrie, 3.0: 30-7-2013        *
!* R. Bannister, Oct 2013           *
!* R. Bannister, vn1.4da 20-10-2017 *
!*                                  *
!************************************

USE DefConsTypes, ONLY :   &
    UM_type,               &
    dims_type,             &
    ABC_type,              &
    ZREAL8,                &
    nlongs, nlevs,         &
    C, f,                  &
    Re, dx, recipdx,       &
    recip2dx, g, theta_r,  &
    Rd, kappa, p00

IMPLICIT NONE

INCLUDE "Boundaries.interface"


! Declare Top level variables
TYPE(UM_type),     INTENT(IN)         :: um_data
TYPE(ABC_type),    INTENT(INOUT)      :: state
TYPE(dims_type),   INTENT(IN)         :: dims

! Declare local variables
REAL(ZREAL8)                          :: cons, Re2
REAL(ZREAL8)                          :: rhopsumtemp, rhopsum(1:nlevs)
REAL(ZREAL8)                          :: rhobound(1:nlevs)
REAL(ZREAL8)                          :: diff(1:nlevs+1), ux, uxm1, px, px1, pxm1, pxp1, deltaz, dpdz
REAL(ZREAL8)                          :: theta_sum, theta_z(1:nlevs)
REAL(ZREAL8)                          :: pprimed(1:nlongs,1:nlevs), pav
REAL(ZREAL8)                          :: rhoprimed(1:nlongs,1:nlevs), rhoav(1:nlevs), rho0
REAL(ZREAL8)                          :: C_um
REAL(ZREAL8)                          :: ratio, ratioav, stddev
REAL(ZREAL8)                          :: INT_HF
INTEGER                               :: x, z, i, j, k


! Estimate C from the UM (method 1 - the formula given in the paper Sec 2.1.2, point 4)
! -------------------------------------------------------------------------------------
C_um = Rd * theta_r * SUM(um_data % exner_pressure(1:nlongs,1:nlevs)) / &
        ( (1.-kappa) * REAL(nlongs*nlevs) )
PRINT *, 'Value of C found from UM exner data is ', C_um

! Estimate C from the UM (method 2 - using p'/rho')
! -------------------------------------------------

! Calculate pprimed
pprimed(1:nlongs,1:nlevs) = p00 * (um_data % exner_pressure(1:nlongs,1:nlevs))**(1./kappa)
DO z = 1, nlevs
  pav                     = SUM(pprimed(1:nlongs,z)) / REAL(nlongs)
  pprimed(1:nlongs,z)     = pprimed(1:nlongs,z) - pav
END DO

! Calculate rhoprimed
rhoprimed(1:nlongs,1:nlevs) = um_data % density(1:nlongs,1:nlevs) / (Re * Re)
DO z = 1, nlevs
  rhoav(z)                  = SUM(rhoprimed(1:nlongs,z)) / REAL(nlongs)
  rhoprimed(1:nlongs,z)     = rhoprimed(1:nlongs,z) - rhoav(z)
END DO

! Calculate average ratio p'/rho'
ratioav = 0.
DO x = 1, nlongs
  DO z = 1, nlevs
    ratio   = pprimed(x,z) / rhoprimed(x,z)
    ratioav = ratioav + ABS(ratio)
  END DO
END DO
ratioav = ratioav / REAL(nlongs * nlevs)

! Calculate standard deviation of this
stddev = 0
DO x = 1, nlongs
  DO z = 1, nlevs
    ratio  = pprimed(x,z) / rhoprimed(x,z) - ratioav
    stddev = stddev + ratio * ratio
  END DO
END DO
stddev = SQRT(stddev / REAL(nlongs * nlevs))
PRINT *, 'Value of C found from p-primed/rho-primed is ', ratioav
PRINT *, '  with standard deviation ', stddev


! Caculate global average rho
rho0 = SUM(rhoav(1:nlevs)) / REAL(nlevs)
PRINT *, 'Global average rho (for rho0) = ', rho0

! Set u from UM
!--------------
PRINT*, 'Setting u ...'
!state % u(1:nlongs, 1:nlevs) = state % u(1:nlongs, 1:nlevs) + um_data % u(1:nlongs, 1:nlevs)
state % u(1:nlongs, 1:nlevs) = um_data % u(1:nlongs, 1:nlevs)
CALL BoundaryMod (state % u(1:nlongs, 1:nlevs))
DO z =1, nlevs
  cons                  = SUM(state % u(1:nlongs, z)) / REAL(nlongs)
  state % u(1:nlongs,z) = state % u(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_u=.TRUE.)

! Set v from UM and enforce integral of v dx to be 0 on each vert level
! ----------------------------------------------------------------------
PRINT*, 'Setting v ...'
!state % v(1:nlongs, 1:nlevs) = state % v(1:nlongs, 1:nlevs) + um_data % v(1:nlongs, 1:nlevs)
state % v(1:nlongs, 1:nlevs) = um_data % v(1:nlongs, 1:nlevs)
CALL BoundaryMod (state % v(1:nlongs, 1:nlevs))
DO z =1, nlevs
  cons                  = SUM(state % v(1:nlongs, z)) / REAL(nlongs)
  state % v(1:nlongs,z) = state % v(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_v=.TRUE.)


! Set p by integrating the geostrophic balance eq and enforce integral of p dx to be 0 on each vert level
! -------------------------------------------------------------------------------------------------------
PRINT*, 'Setting p by integrating geostrophic balance ...'
DO z = 1, nlevs
  !state % r(1,z) = 0.
  DO x = 2, nlongs
    state % r(x,z) = state % r(x-1,z) + f * dx * (state % v(x,z) + state % v(x-1,z)) / (2. * C)
  END DO
  cons                  = SUM(state % r(1:nlongs,z)) / REAL(nlongs)
  state % r(1:nlongs,z) = state % r(1:nlongs,z) - cons
END DO
CALL Boundaries (state, set_r=.TRUE.)

! Set rho
! -------
PRINT*, 'Setting rho ...'
state % rho(1:nlongs,1:nlevs) = state % r(1:nlongs,1:nlevs) + 1.0
CALL Boundaries (state, set_rho=.TRUE.)

! Set bp from the hydrostatic relationship
! ----------------------------------------
PRINT*, 'Setting bp from hydrostatic balance ...'
DO z = 1, nlevs
  DO x = 1, nlongs
    state % b(x,z) = C * (state % r(x,z+1) - state % r(x,z)) /      &
                         (Dims % half_levs(z+1) - Dims % half_levs(z))
  END DO
END DO
DO z =1, nlevs
  cons                  = SUM(state % b(1:nlongs, z)) / REAL(nlongs)
  state % b(1:nlongs,z) = state % b(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_b=.TRUE.)

! Calculate w from the continuity equation
! ----------------------------------------
PRINT*, 'Setting w from the continuity equation ...'
DO x = 1, nlongs
  ! Calculate the d(rho u)/dx for a column at this longitude
  DO z = 1, nlevs
    ! u interpolated to full-level
    ux      = INT_HF(state % u(x,z), state % u(x,z+1), z, Dims)
    uxm1    = INT_HF(state % u(x-1,z), state % u(x-1,z+1), z, Dims)
    ! rho interpolated to full-level
    px    = INT_HF(state % rho(x,z), state % rho(x,z+1), z, Dims)
    pxm1  = INT_HF(state % rho(x-1,z), state % rho(x-1,z+1), z, Dims)
    pxp1  = INT_HF(state % rho(x+1,z), state % rho(x+1,z+1), z, Dims)
    ! d(p u)/dx
    diff(z) = recip2dx * ( (pxp1 + px) * ux - (px + pxm1) * uxm1 )
  END DO

  ! Linear extrapolation of diff to nlevs + 1
  diff(nlevs+1) = diff(nlevs) +                                              &
                  ( diff(nlevs) - diff(nlevs-1) ) *                          &
                  ( Dims % full_levs(nlevs+1) - Dims % full_levs(nlevs) ) /  &
                  ( Dims % full_levs(nlevs) - Dims % full_levs(nlevs-1) )

  ! Integrate the continuity equation from the ground upwards
  !state % w(x,0) = 0.


  ! Following loops are alternative codes (one should be commented out)
  ! Cannot remember justification for one over the other

!  ! First consider state % w as rho*w
!  DO z = 1, nlevs - 1
!   ! Level width
!    deltaz  = Dims % half_levs(z+1) - Dims % half_levs(z)
!    ! Integrate rho*w
!    state % w(x,z) = state % w(x,z-1) - diff(z) * deltaz
!  END DO
!  !state % w(x,nlevs) = 0.0
!  ! Second, convert rho*w to just w
!  DO z = 1, nlevs - 1
!    px             = INT_HF(state % rho(x,z), state % rho(x,z+1), z, Dims)
!    state % w(x,z) = state % w(x,z) / px
!  END DO

  DO z = 1, nlevs - 1
    ! rho interpolated to full-level
    px = INT_HF(state % rho(x,z), state % rho(x,z+1), z, Dims)
    ! Level width
    deltaz  = Dims % half_levs(z+1) - Dims % half_levs(z)
    ! dp/dz
    dpdz  = (state % rho(x,z+1) - state % rho(x,z)) / deltaz
    ! Increment w
    state % w(x,z) = 0.0 * state % w(x,z-1) - (diff(z) + dpdz * state % w(x,z-1)) * deltaz / px
  END DO

END DO

DO z =1, nlevs
  cons                  = SUM(state % w(1:nlongs, z)) / REAL(nlongs)
  state % w(1:nlongs,z) = state % w(1:nlongs, z) - cons
END DO
CALL Boundaries (state, set_w=.TRUE.)



END SUBROUTINE Process_um_data
