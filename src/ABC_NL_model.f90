SUBROUTINE ABC_NL_model (statein, stateout, Dims)

! Forward model integrated over dt
! adjust is the input state

USE DefConsTypes, ONLY :               &
  ZREAL8,                              &
  ABC_type,                            &
  dims_type,                           &
  Averages_type,                       &
  nlongs, nlevs,                       &
  bdiva_f, deltat, recip_alpha_f,      &
  recipdx, f, recip2dx, bdiva_N, half, &
  recip_alpha_N, A, B, C, dt,          &
  Adv_tracer

IMPLICIT NONE

INCLUDE "Boundaries.interface"

! Declare global parameters
!--------------------------
TYPE(ABC_type),  INTENT(IN)    :: statein
TYPE(ABC_type),  INTENT(INOUT) :: stateout
TYPE(dims_type), INTENT(IN)    :: Dims

! Declare local variables
!------------------------
TYPE(ABC_type)                 :: adjust, adjust_p1
TYPE(Averages_type)            :: aves
REAL(ZREAL8)                   :: u_at_b(1:nlongs,1:nlevs), w_at_u(1:nlongs,1:nlevs)
REAL(ZREAL8)                   :: u_at_v(1:nlongs,1:nlevs), w_at_v(1:nlongs,1:nlevs)
REAL(ZREAL8)                   :: upstr_u, upstr_w, rhodivu
REAL(ZREAL8)                   :: deltat_recip_alpha_f, half_f, half_deltat_C
REAL(ZREAL8)                   :: deltat_recip_alpha_N, deltat_recip_alpha_N_N2, deltat_B
REAL(ZREAL8)                   :: deltat_recip_alpha_f_f, B_dt
REAL(ZREAL8)                   :: drdx1, drdx2, drdz, Nsq

! Declare functions
!------------------
REAL(ZREAL8)                  :: INT_FH, INT_HF

! Declare Flags
!--------------
LOGICAL                       :: flag_u, flag_w, flag_u_at_b, flag_w_at_b
LOGICAL                       :: flag_w_at_u, flag_u_at_v, flag_w_at_v

! Declare Loop counters
!----------------------
INTEGER                       :: i, k, tsmall

! Initialise
CALL Initialise_Averages (aves)
Nsq                     = A * A
deltat_recip_alpha_f    = deltat * recip_alpha_f
half_f                  = half * f
half_deltat_C           = half * deltat * C
deltat_recip_alpha_N    = deltat * recip_alpha_N
deltat_recip_alpha_N_N2 = deltat * recip_alpha_N * Nsq
deltat_B                = deltat * B
deltat_recip_alpha_f_f  = deltat * recip_alpha_f * f
B_dt                    = B * dt

! copy input state into adjust array
adjust = statein

!****************************
!***** ADJUSTMENT STAGE *****
!****************************

!Begin temporal loop, with small timestep
DO tsmall = 1, 2

!*** Forward-backward scheme ***

  DO k = 1, nlevs
    DO i = 1, nlongs
      drdx1              = ( adjust % r(i+1,k) - adjust % r(i,k) ) * recipdx
      drdx2              = ( adjust % r(i+1,k) - adjust % r(i-1,k) ) * recip2dx

      adjust_p1 % u(i,k) = bdiva_f * adjust % u(i,k) +                                                &
                           deltat_recip_alpha_f * ( half_f * (adjust % v(i,k) + adjust % v(i+1,k) ) - &
                                                    C * drdx1 )

      adjust_p1 % v(i,k) = bdiva_f * adjust % v(i,k) +                                                &
                           deltat_recip_alpha_f_f * ( half_deltat_C * drdx2 -                         &
                                                      half * ( adjust % u(i-1,k) + adjust % u(i,k) ) )

      drdz               = ( adjust % r(i,k+1) - adjust % r(i,k) ) * dims % recip_half_kp1_k(k)

      adjust_p1 % w(i,k) = bdiva_N * adjust % w(i,k) +                                                &
                           deltat_recip_alpha_N * ( adjust % b(i,k) - C * drdz )

      adjust_p1 % b(i,k) = bdiva_N * adjust % b(i,k) +                                                &
                           deltat_recip_alpha_N_N2 * ( half_deltat_C * drdz - adjust % w(i,k) )
    END DO
  END DO

  ! Reassign variables
  !-------------------
  adjust % u(1:nlongs,1:nlevs) = adjust_p1 % u(1:nlongs,1:nlevs)
  adjust % v(1:nlongs,1:nlevs) = adjust_p1 % v(1:nlongs,1:nlevs)
  adjust % w(1:nlongs,1:nlevs) = adjust_p1 % w(1:nlongs,1:nlevs)
  adjust % b(1:nlongs,1:nlevs) = adjust_p1 % b(1:nlongs,1:nlevs)
  adjust % w(1:nlongs,0)       = adjust_p1 % w(1:nlongs,0)

  CALL Boundaries (adjust, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_b=.TRUE.)

  !Store all adjustment variables
  !------------------------------
  IF (tsmall == 1) THEN
    aves % u_1(0:nlongs+1,0:nlevs+1) = adjust % u(0:nlongs+1,0:nlevs+1)
    aves % w_1(0:nlongs+1,0:nlevs+1) = adjust % w(0:nlongs+1,0:nlevs+1)
  ELSE IF (tsmall == 2 ) THEN
    aves % u_2(0:nlongs+1,0:nlevs+1) = adjust % u(0:nlongs+1,0:nlevs+1)
    aves % w_2(0:nlongs+1,0:nlevs+1) = adjust % w(0:nlongs+1,0:nlevs+1)
  END IF


!!!!!  adjust % rho(1:nlongs,1:nlevs) = 1.0 + adjust % r(1:nlongs,1:nlevs)
!!!!!
!!!!!  !Apply boundary conditions to full rho
!!!!!  CALL Boundaries (adjust, set_rho=.TRUE.)

  !Calculate rnext using a forward upstream
  !----------------------------------------

  !Calculate advecting interpolations
  DO k = 1, nlevs
    DO i = 1, nlongs
      u_at_v(i,k) = half * ( adjust % u(i-1, k) + adjust % u(i,k) )
      w_at_v(i,k) = INT_FH( adjust % w(i,k-1), adjust % w(i,k), k , Dims )
    END DO
  END DO

  DO k = 1, nlevs
    DO i = 1,nlongs

      !Calculate flags
      flag_u_at_v = ( u_at_v(i,k) .LE. 0.0 )
      flag_w_at_v = ( w_at_v(i,k) .LE. 0.0 )

      !Calculate upstream terms
      IF (flag_u_at_v) THEN
        upstr_u = u_at_v(i,k) * ( adjust % rho(i+1,k) - adjust % rho(i,k) ) * recipdx
      ELSE
        upstr_u = u_at_v(i,k) * ( adjust % rho(i,k) - adjust % rho(i-1,k) ) * recipdx
      END IF

      IF (flag_w_at_v) THEN
        upstr_w = w_at_v(i,k) * ( adjust % rho(i,k+1) - adjust % rho(i,k) ) *        &
                                Dims % recip_half_kp1_k(k)
      ELSE
        upstr_w = w_at_v(i,k) * ( adjust % rho(i,k) - adjust % rho(i,k-1) ) *        &
                                Dims % recip_half_k_km1(k)
      END IF

      rhodivu  =  adjust % rho(i,k) *                                                &
                  ( adjust % u(i,k) - adjust % u(i-1,k) ) * recipdx +                &
                  ( adjust % w(i,k) - adjust % w(i,k-1) ) * Dims % recip_full_k_km1(k)

      adjust_p1 % r(i,k) = adjust % r(i,k) - deltat_B *    &
                           ( rhodivu +                     &   ! rho div u
                             upstr_w + upstr_u )               ! u dot grad rho
    END DO
  END DO

  ! Update fields
  adjust % r(1:nlongs,1:nlevs)   = adjust_p1 % r(1:nlongs,1:nlevs)
  adjust % rho(1:nlongs,1:nlevs) = 1.0 + adjust % r(1:nlongs,1:nlevs)
  ! Apply boundary conditions
  CALL Boundaries (adjust, set_r=.TRUE., set_rho=.TRUE.)

  !End of temporal loop for tsmall
END DO

!*************************************
!***** Calculate advecting means *****
!*************************************

aves % u_m(1:nlongs,1:nlevs) = 0.5 * ( aves % u_1(1:nlongs,1:nlevs) + aves % u_2(1:nlongs,1:nlevs) )
aves % w_m(1:nlongs,1:nlevs) = 0.5 * ( aves % w_1(1:nlongs,1:nlevs) + aves % w_2(1:nlongs,1:nlevs) )

! Apply Boundary Conditions to advecting means
!---------------------------------------------
aves % u_m(0,0:nlevs+1)         = aves % u_m(nlongs, 0:nlevs+1)
aves % w_m(0,0:nlevs+1)         = aves % w_m(nlongs, 0:nlevs+1)
aves % u_m(nlongs+1,0:nlevs+1)  = aves % u_m(1, 0:nlevs+1)
aves % w_m(nlongs+1,0:nlevs+1)  = aves % w_m(1, 0:nlevs+1)
aves % u_m( 0:nlongs+1,nlevs+1) = aves % u_m( 0:nlongs+1, nlevs)
aves % u_m( 0:nlongs+1,0)       = -1.0 * aves % u_m( 0:nlongs+1,1)
aves % w_m(0:nlongs+1,0)        = 0.0
aves % w_m(0:nlongs+1,nlevs)    = 0.0
aves % w_m(0:nlongs+1,nlevs+1)  = 0.0

!Reassign r
stateout % r(1:nlongs,1:nlevs) = adjust % r(1:nlongs,1:nlevs)
CALL Boundaries (stateout, set_r=.TRUE.)


!***************************
!***** ADVECTION STAGE *****
!***************************

!Calculate advecting interpolations
DO k = 1, nlevs
  DO i = 1, nlongs
    u_at_b(i,k) = half * ( INT_HF(aves % u_m(i-1,k), aves % u_m(i-1,k+1), k, Dims)  + &
                           INT_HF(aves % u_m(i,k), aves % u_m(i,k+1), k, Dims) )
    u_at_v(i,k) = half * ( aves % u_m(i-1, k) + aves % u_m(i,k) )
    w_at_u(i,k) = half * ( INT_FH( aves % w_m(i,k-1), aves % w_m(i, k), k , Dims)   + &
                           INT_FH( aves % w_m(i+1,k-1), aves % w_m(i+1, k), k , Dims))
    w_at_v(i,k) = INT_FH( aves % w_m(i,k-1), aves % w_m(i, k), k , Dims)
  END DO
END DO

DO k = 1, nlevs
  DO i = 1,nlongs

    !Calculate flags
    !---------------
    flag_u      = ( aves % u_m(i,k) .LE. 0.0 )
    flag_w      = ( aves % w_m(i,k) .LE. 0.0 )
    flag_u_at_b = ( u_at_b(i,k) .LE. 0.0 )
    flag_w_at_u = ( w_at_u(i,k) .LE. 0.0 )
    flag_u_at_v = ( u_at_v(i,k) .LE. 0.0 )
    flag_w_at_v = ( w_at_v(i,k) .LE. 0.0 )

    ! u- component
    !-------------
    IF (flag_u) THEN
      upstr_u = aves % u_m(i,k) * ( statein % u(i+1,k) - statein % u(i,k) ) * recipdx
    ELSE
      upstr_u = aves % u_m(i,k) * ( statein % u(i,k) - statein % u(i-1,k) ) * recipdx
    END IF
    IF (flag_w_at_u) THEN
      upstr_w = w_at_u(i,k) * ( statein % u(i,k+1) - statein % u(i,k) ) *                 &
                              Dims % recip_half_kp1_k(k)
    ELSE
      upstr_w = w_at_u(i,k) * ( statein % u(i,k) - statein % u(i,k-1) ) *                 &
                              Dims % recip_half_k_km1(k)
    END IF
    stateout % u(i,k) = adjust % u(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad u

    ! v- component
    !-------------
    IF (flag_u_at_v) THEN
      upstr_u = u_at_v(i,k) * ( statein % v(i+1,k) - statein % v(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_v(i,k) * ( statein % v(i,k) - statein % v(i-1,k) ) * recipdx
    END IF
    IF (flag_w_at_v) THEN
      upstr_w = w_at_v(i,k) * ( statein % v(i,k+1) - statein % v(i,k) ) *                &
                              Dims % recip_half_kp1_k(k)
    ELSE
      upstr_w = w_at_v(i,k) * ( statein % v(i,k) - statein % v(i,k-1) ) *                &
                              Dims % recip_half_k_km1(k)
    END IF

    stateout % v(i,k) = adjust % v(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad v

    ! w- component
    !-------------
    IF (flag_u_at_b) THEN
      upstr_u = u_at_b(i,k) * ( statein % w(i+1,k) - statein % w(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_b(i,k) * ( statein % w(i,k) - statein % w(i-1,k) ) * recipdx
    END IF
    IF (flag_w) THEN
      upstr_w = aves % w_m(i,k) *                                                        &
                ( statein % w(i,k+1) - statein % w(i,k) ) *                              &
                Dims % recip_full_kp1_k(k)
    ELSE
      upstr_w = aves % w_m(i,k) *                                                        &
                ( statein % w(i,k) - statein % w(i,k-1) ) *                              &
                Dims % recip_full_k_km1(k)
    END IF

    stateout % w(i,k) = adjust % w(i,k) - B_dt * ( upstr_u + upstr_w )   ! vec(u) . grad w


    ! b' component
    !-------------
    IF (flag_u_at_b) THEN
      upstr_u = u_at_b(i,k) * ( statein % b(i+1,k) - statein % b(i,k) ) * recipdx
    ELSE
      upstr_u = u_at_b(i,k) * ( statein % b(i,k) -  statein % b(i-1,k) ) * recipdx
    END IF
    IF (flag_w) THEN
      upstr_w = aves % w_m(i,k) *                                                        &
                ( statein % b(i,k+1) - statein % b(i,k) ) *                              &
                Dims % recip_full_kp1_k(k)
    ELSE
      upstr_w = aves % w_m(i,k) *                                                        &
                ( statein % b(i,k) - statein % b(i,k-1) ) *                              &
                Dims % recip_full_k_km1(k)

    END IF

    stateout % b(i,k) = adjust % b(i,k) - B_dt * ( upstr_u + upstr_w )   ! u . grad b
  END DO
END DO

! Apply boundary conditions, u, v, w, b'
CALL Boundaries (stateout, set_u=.TRUE., set_v=.TRUE., set_w=.TRUE., set_b=.TRUE.)


!*************************
!***** Advect Tracer *****
!*************************

IF (Adv_tracer) THEN
  DO k = 1, nlevs
    DO i = 1,nlongs

      !Calculate flags
      flag_u_at_b = ( u_at_b(i,k) .LE. 0.0 )
      flag_w_at_b = ( aves % w_m(i,k) .LE. 0.0 )

      IF (flag_u_at_b) THEN
        upstr_u = u_at_b(i,k) * ( statein % tracer(i+1,k) - statein % tracer(i,k) ) * recipdx
      ELSE
        upstr_u = u_at_b(i,k) * ( statein % tracer(i,k) - statein % tracer(i-1,k) ) * recipdx
      END IF
      IF (flag_w_at_b) THEN
        upstr_w = aves % w_m(i,k) * ( statein % tracer(i,k+1) - statein % tracer(i,k) ) *    &
                  Dims % recip_full_kp1_k(k)
      ELSE
        upstr_w = aves % w_m(i,k) * ( statein % tracer(i,k) - statein % tracer(i,k-1) ) *    &
                  Dims % recip_full_k_km1(k)
      END IF

      stateout % tracer (i,k) = statein % tracer(i,k) - dt * ( upstr_u + upstr_w )   ! vec(u) . grad tracer
    END DO
  END DO
  CALL Boundaries (stateout, set_tracer=.TRUE.)
END IF

! Update rho
stateout % rho(1:nlongs,1:nlevs) = 1.0 + stateout % r(1:nlongs,1:nlevs)
CALL Boundaries (stateout, set_rho=.TRUE.)

END SUBROUTINE ABC_NL_model
