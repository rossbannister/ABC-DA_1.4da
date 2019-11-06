SUBROUTINE ABC_NL_ModelDriver_DA (States, dims, ntimesteps, steps4da)

!********************************************************************************
!*                                                                              *
!*  Driver Routine for Non linear forward model for use in DA                   *
!*                                                                              *
!*  StateA                 - initial condition (modified on output)             *
!*  dims                   - dimension data                                     *
!*  ntimesteps             - total number of standard model timesteps           *
!*  steps4da               - number of timesteps per DA timestep within the DA  *
!*                           window                                             *
!*                                                                              *
!*   R. Bannister, 1.4da 25-01-2018                                             *
!*                                                                              *
!********************************************************************************

USE DefConsTypes, ONLY :         &
  ZREAL8,                        &
  ABC_type,                      &
  dims_type,                     &
  nlongs,                        &
  nlevs,                         &
  dt, dt_da, deltat, dx, dz,     &
  Lengthscale_diagnostics,&
  A, B, C, f, g

IMPLICIT NONE

INCLUDE "Boundaries.interface"

! Parameters
!-----------
TYPE(ABC_type),  INTENT(INOUT) :: States(0:steps4da)
TYPE(dims_type), INTENT(IN)    :: dims
INTEGER,         INTENT(IN)    :: ntimesteps
INTEGER,         INTENT(IN)    :: steps4da

! Local Variables
!----------------
INTEGER                        :: t, da_timestep, stepsperda
TYPE(ABC_type)                 :: StateA, StateB


IF ( (dt/deltat /= 2.0 )) THEN
  PRINT*, '******************************'
  PRINT*, '********** Error *************'
  PRINT*, '  dt /deltat .NE. 2'
  STOP
ENDIF

stepsperda = ntimesteps / steps4da

PRINT*, 'Model timestep         ', dt
PRINT*, 'DA timestep            ', dt_da
PRINT*, 'Total no. of timesteps ', ntimesteps
PRINT*, 'Total run length       ', dt*ntimesteps,' s'
PRINT*, 'Total run length       ', dt*ntimesteps/3600,' hrs'
PRINT*, 'Number of DA time steps', steps4da

! Make sure that the rho field is consistent with the r field
States(0) % rho(1:nlongs,1:nlevs) = 1.0 + States(0) % r(1:nlongs,1:nlevs)

! Apply boundary conditions
CALL Boundaries (States(0), set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_r=.TRUE., &
                            set_b=.TRUE., set_rho=.TRUE., set_tracer=.TRUE.)


PRINT*,'----------------------------'
PRINT*,'    Running prognostic model'
PRINT*,'----------------------------'

StateA      = States(0)
da_timestep = 0


DO t = 1, ntimesteps
  !PRINT*,' > ', t, ' of ', ntimesteps


  CALL ABC_NL_model (StateA, StateB, dims)


  ! At da times, store state
  !-------------------------
  IF ( (stepsperda * INT(REAL(t)/REAL(stepsperda))) == t) THEN
    da_timestep         = da_timestep + 1
    ! PRINT *, 'Storing this timestep for DA', da_timestep
    States(da_timestep) = StateB
  END IF

  ! Re-assign
  StateA = StateB

END DO



END SUBROUTINE ABC_NL_ModelDriver_DA
