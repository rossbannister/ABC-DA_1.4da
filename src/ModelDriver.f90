SUBROUTINE ABC_NL_ModelDriver ( StateIn, dims, ntimesteps, ndumps,   &
                                fileout, diagnostics_file )

!********************************************************************************
!*                                                                              *
!*  Driver Routine for Non linear forward model                                 *
!*                                                                              *
!*  StateIn                - output state of ABC_type                           *
!*  dims                   - dimension data                                     *
!*  ntimesteps             - total number of timesteps for model integration    *
!*  ndumps                 - number of output times                             *
!*  fileout                - output file                                        *
!*  diagnostics_file       - diagnostics file                                   *
!*                                                                              *
!*   R. Petrie, 2.0:  10-6-2011                                                 *
!*   R. Petrie, 3.0:  30-7-2013                                                 *
!*   R. Bannister, 1.4da 22-10-2017                                             *
!*                                                                              *
!********************************************************************************

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  dims_type,              &
  nlongs,                 &
  nlevs,                  &
  dt, deltat, dx, dz,     &
  Lengthscale_diagnostics,&
  A, B, C, f, g

IMPLICIT NONE

INCLUDE "Boundaries.interface"

! Parameters
!-----------
TYPE(ABC_type),  INTENT(IN) :: StateIn
TYPE(dims_type), INTENT(IN) :: dims
INTEGER,         INTENT(IN) :: ntimesteps, ndumps
CHARACTER(*),    INTENT(IN) :: fileout
CHARACTER(*),    INTENT(IN) :: diagnostics_file

! Local Variables
!----------------
INTEGER                     :: dumpspacing
INTEGER                     :: t, diagnostics_unit
REAL(ZREAL8)                :: L_horiz_u, L_horiz_v, L_horiz_w, L_horiz_r, L_horiz_b, L_horiz_uv, L_horiz_rho
REAL(ZREAL8)                :: L_vert_u, L_vert_v, L_vert_w, L_vert_r, L_vert_b, L_vert_uv, L_vert_rho
REAL(ZREAL8)                :: rms_u, rms_v, rms_w, rms_r, rms_b, rms_uv, rms_rho, au
REAL(ZREAL8)                :: ca, cg, Ro, M, Fr, RrA, Rrg, ar, sr, rr

TYPE(ABC_type),  INTENT(IN) :: StateA, StateB

! Functions
! ---------
REAL(ZREAL8)                :: RMS


! Consistency checks
dumpspacing = ntimesteps / ndumps

IF ( (dt/deltat /= 2.0 )) THEN
  PRINT*, '******************************'
  PRINT*, '********** Error *************'
  PRINT*, '  dt /deltat .NE. 2'
  STOP
ENDIF

PRINT*, 'Total no. of timesteps ', ntimesteps
PRINT*, 'Total run length       ', dt*ntimesteps,' s'
PRINT*, 'Total run length       ', dt*ntimesteps/3600,' hrs'
PRINT*, 'Number of output times ', ndumps
PRINT*, 'Output times every     ', dt*ntimesteps/(3600*ndumps),' hrs'


! Diagnostics output file
OPEN (diagnostics_unit, file=diagnostics_file)

StateA = StateIn

! Apply boundary conditions
CALL Boundaries (StateA, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_r=.TRUE., &
                         set_b=.TRUE., set_rho=.TRUE., set_tracer=.TRUE.)

! Calculate some diagnostics and write Initial Conditions
!-------------------------
CALL Calc_geost(StateA, dims)
CALL Calc_hydro(StateA, dims)
CALL Calc_vert_mom_source(StateA, dims)
CALL Calc_horiz_div(StateA, dims)
CALL Calc_horiz_vort(StateA, dims)
CALL Energy(StateA)
CALL Effective_buoyancy(StateA, dims)
CALL Boundaries (StateA, set_beff=.TRUE.)
CALL Write_state_2d (fileout, StateA, dims, ndumps+1, 0, dumpspacing)



PRINT*,'----------------------------'
PRINT*,'    Running prognostic model'
PRINT*,'----------------------------'

DO t = 1, ntimesteps
 ! PRINT*,' > ', t, ' of ', ntimesteps


  CALL ABC_NL_model (StateA, StateB, dims)


  ! Calculate the energy contributions
  CALL Energy (StateB)

  ! At output times
  !-----------------
  IF ((t/dumpspacing)*dumpspacing == t) THEN
    ! Calculate diagnostics
    !------------------------
    CALL Calc_hydro (StateB, dims)
    CALL Calc_geost (StateB, dims)
    CALL Calc_vert_mom_source(StateB, dims)
    CALL Calc_horiz_div(StateB, dims)
    CALL Calc_horiz_vort(StateB, dims)
    CALL Effective_buoyancy (StateB, dims)
    CALL Write_state_2d (fileout, StateB, dims, ndumps+1, t/dumpspacing, dumpspacing)
  ENDIF

  WRITE (diagnostics_unit, *) t, t*dt, StateB % Kinetic_Energy,      &
                                       StateB % Buoyant_Energy,      &
                                       StateB % Elastic_Energy,      &
                                       StateB % Total_Energy

  ! Lengthscale diagnostics at the last timestep
  IF (Lengthscale_diagnostics .AND. (t == ntimesteps)) THEN
    WRITE (diagnostics_unit, *) '# -----------------------------------------'
    WRITE (diagnostics_unit, *) '# Diagnostics for last step (SI units)'
    WRITE (diagnostics_unit, *) '# -----------------------------------------'

    WRITE (*,*) 'Diagnosing lengthscales for u ...'
    CALL Lscales_from_fft(StateB % u(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_u, L_vert_u)
    L_horiz_u = L_horiz_u * dx
    L_vert_u  = L_vert_u * dz
    WRITE (*,*) 'Diagnosing RMS for u ...'
    CALL Magnitude_rms(StateB % u(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_u)

    WRITE (*,*) 'Diagnosing lengthscales for v ...'
    CALL Lscales_from_fft(StateB % v(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_v, L_vert_v)
    L_horiz_v = L_horiz_v * dx
    L_vert_v  = L_vert_v * dz
    WRITE (*,*) 'Diagnosing RMS for v ...'
    CALL Magnitude_rms(StateB % v(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_v)

    WRITE (*,*) 'Diagnosing lengthscales for w ...'
    CALL Lscales_from_fft(StateB % w(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_w, L_vert_w)
    L_horiz_w = L_horiz_w * dx
    L_vert_w  = L_vert_w * dz
    WRITE (*,*) 'Diagnosing RMS for w ...'
    CALL Magnitude_rms(StateB % w(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_w)

    WRITE (*,*) 'Diagnosing lengthscales for r ...'
    CALL Lscales_from_fft(StateB % r(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_r, L_vert_r)
    L_horiz_r = L_horiz_r * dx
    L_vert_r  = L_vert_r * dz
    WRITE (*,*) 'Diagnosing RMS for r ...'
    CALL Magnitude_rms(StateB % r(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_r)

    WRITE (*,*) 'Diagnosing lengthscales for rho ...'
    CALL Lscales_from_fft(StateB % rho(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_rho, L_vert_rho)
    L_horiz_rho = L_horiz_rho * dx
    L_vert_rho  = L_vert_rho * dz
    WRITE (*,*) 'Diagnosing RMS for rho ...'
    CALL Magnitude_rms(StateB % rho(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_rho)

    WRITE (*,*) 'Diagnosing lengthscales for b ...'
    CALL Lscales_from_fft(StateB % b(1:nlongs,2:nlevs-1), nlongs, nlevs-2, L_horiz_b, L_vert_b)
    L_horiz_b = L_horiz_b * dx
    L_vert_b  = L_vert_b * dz
    WRITE (*,*) 'Diagnosing RMS for b ...'
    CALL Magnitude_rms(StateB % b(1:nlongs,2:nlevs-1), nlongs, nlevs-2, rms_b)

    ! Average data for horizontal wind
    L_horiz_uv = (L_horiz_u + L_horiz_v) / 2.0
    L_vert_uv  = (L_vert_u + L_vert_v) / 2.0
    rms_uv     = (rms_u + rms_v) / 2.0


    ! Acoustic wave speed
    ca  = SQRT(B*C)
    ! Gravity wave speed
    cg  = SQRT(g*L_vert_r)
    ! Rossby number
    Ro  = rms_u / (f * L_horiz_r)
    ! Mach number
    M   = rms_u / ca
    ! Froude number
    Fr  = rms_u / cg
    ! Rossby radius calculated with static stability
    RrA = A * L_vert_r / f
    ! Rossby radius calculated with gravity wave speed
    Rrg = cg / f
    ! Aspect ratio
    ar  = L_vert_u / L_horiz_u
    ! Speed ratio
    sr  = rms_w / rms_u
    ! Ratio of ratios
    rr  = ar / sr
    ! Ratio of horizontal wind speeds
    au = rms_v / rms_u

    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'A,                         ', A
    WRITE (diagnostics_unit, *) 'B,                         ', B
    WRITE (diagnostics_unit, *) 'C,                         ', C
    WRITE (diagnostics_unit, *) 'f,                         ', f
    WRITE (diagnostics_unit, *) 'g,                         ', g
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'u L_h,                     ', L_horiz_u
    WRITE (diagnostics_unit, *) 'u L_v,                     ', L_vert_u
    WRITE (diagnostics_unit, *) 'u RMS,                     ', rms_u
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'v L_h,                     ', L_horiz_v
    WRITE (diagnostics_unit, *) 'v L_v,                     ', L_vert_v
    WRITE (diagnostics_unit, *) 'v RMS,                     ', rms_v
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'hor-sp L_h,                ', L_horiz_uv
    WRITE (diagnostics_unit, *) 'hor-sp L_v,                ', L_vert_uv
    WRITE (diagnostics_unit, *) 'hor-sp RMS,                ', rms_uv
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'w L_h,                     ', L_horiz_w
    WRITE (diagnostics_unit, *) 'w L_v,                     ', L_vert_w
    WRITE (diagnostics_unit, *) 'w RMS,                     ', rms_w
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'r L_h,                     ', L_horiz_r
    WRITE (diagnostics_unit, *) 'r L_v,                     ', L_vert_r
    WRITE (diagnostics_unit, *) 'r RMS,                     ', rms_r
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'rho L_h,                   ', L_horiz_rho
    WRITE (diagnostics_unit, *) 'rho L_v,                   ', L_vert_rho
    WRITE (diagnostics_unit, *) 'rho RMS,                   ', rms_rho
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'b L_h,                     ', L_horiz_b
    WRITE (diagnostics_unit, *) 'b L_v,                     ', L_vert_b
    WRITE (diagnostics_unit, *) 'b RMS,                     ', rms_b
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'ca,                        ', ca
    WRITE (diagnostics_unit, *) 'cg,                        ', cg
    WRITE (diagnostics_unit, *) 'RossbyNo,                  ', Ro
    WRITE (diagnostics_unit, *) 'MachNo,                    ', M
    WRITE (diagnostics_unit, *) 'FroudeNo,                  ', Fr
    WRITE (diagnostics_unit, *) 'RbyRad (A),                ', RrA
    WRITE (diagnostics_unit, *) 'RbyRad (cg),               ', Rrg
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'Aspect ratio,              ', ar
    WRITE (diagnostics_unit, *) 'Speed ratio W/U,           ', sr
    WRITE (diagnostics_unit, *) 'Speed ratio V/U            ', au
    WRITE (diagnostics_unit, *) 'aspect ratio / speed ratio ', rr
    WRITE (diagnostics_unit, *) '---'
    WRITE (diagnostics_unit, *) 'MAGNITUDES OF TERMS IN EQUATIONS'
    WRITE (diagnostics_unit, *) 'For a key, see'
    WRITE (diagnostics_unit, *) 'www.met.reading.ac.uk/~ross/DARC/ConvScaleB/ABCmodel/ScaleAnalysis4Paper.html'
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'u1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'u2,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'u3,                        ', B * Ro * sr / ar
    WRITE (diagnostics_unit, *) 'u4,                        ', C * rms_r / (rms_u * f * L_horiz_r)
    WRITE (diagnostics_unit, *) 'u5,                        ', au
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'v1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'v2,                        ', B * Ro * L_horiz_u / L_horiz_v
    WRITE (diagnostics_unit, *) 'v3,                        ', B * Ro * L_horiz_u * sr / L_vert_v
    WRITE (diagnostics_unit, *) 'v4,                        ', 1.0 / au
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'w1,                        ', B * Ro
    WRITE (diagnostics_unit, *) 'w2,                        ', B * Ro * L_horiz_u / L_horiz_w
    WRITE (diagnostics_unit, *) 'w3,                        ', B * Ro * L_horiz_u * sr / L_vert_w
    WRITE (diagnostics_unit, *) 'w4,                        ', C * rms_r / (rms_w * f * L_vert_r)
    WRITE (diagnostics_unit, *) 'w5,                        ', rms_b / (rms_w * f)
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'rhop1,                     ', 1.0
    WRITE (diagnostics_unit, *) 'rho02,                     ', 1.0
    WRITE (diagnostics_unit, *) 'rhop3,                     ', L_horiz_u * sr / L_vert_w
    WRITE (diagnostics_unit, *) ''
    WRITE (diagnostics_unit, *) 'bp1,                       ', B * Ro
    WRITE (diagnostics_unit, *) 'bp2,                       ', B * Ro * L_horiz_u / L_horiz_b
    WRITE (diagnostics_unit, *) 'bp3,                       ', B * Ro * L_horiz_u * sr / L_vert_b
    WRITE (diagnostics_unit, *) 'bp4,                       ', A * A * rms_w / (rms_b * f)
    WRITE (diagnostics_unit, *) '---'

    WRITE (*,*) 'Done diagnostics'
  END IF

  ! Re-assign
  StateA = StateB

END DO

CLOSE (diagnostics_unit)


END SUBROUTINE ABC_NL_ModelDriver
