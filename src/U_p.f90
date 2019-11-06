SUBROUTINE U_p (LS, ControlVar, ModelVar, order,             &
                option_gb, option_hb, option_ab, option_reg, &
                Regression, dims)

! Code to perform the parameter cvt: ModelVar = U_p ControlVar

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  dims_type,              &
  nlongs,                 &
  nlevs

IMPLICIT NONE

INCLUDE "Boundaries.interface"


TYPE(ABC_type),  INTENT(IN)    :: LS
TYPE(CV_type),   INTENT(IN)    :: ControlVar
TYPE(ABC_type),  INTENT(INOUT) :: ModelVar
INTEGER,         INTENT(IN)    :: order
INTEGER,         INTENT(IN)    :: option_gb
INTEGER,         INTENT(IN)    :: option_hb
INTEGER,         INTENT(IN)    :: option_ab
INTEGER,         INTENT(IN)    :: option_reg
REAL(ZREAL8),    INTENT(IN)    :: Regression(1:nlevs, 1:nlevs)
TYPE(dims_type), INTENT(IN)    :: dims

INTEGER                        :: x
REAL(ZREAL8)                   :: r_b(1:nlongs, 1:nlevs)
REAL(ZREAL8)                   :: b_b(1:nlongs, 1:nlevs)
REAL(ZREAL8)                   :: w_b(1:nlongs, 1:nlevs)


IF ((order == 1).OR.(order == 2)) THEN
  ! Traditional kind of control variable transform
  ! ----------------------------------------------

  ! Control variables are psi, chi, (unbalanced) r, (unbalanced) b, (unbalanced) w, tracer

  ! ----------------------------------------------
  ! 1. Compute u and v from psi and chi
  ! ----------------------------------------------
  CALL Helmholtz (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                  ControlVar % v2(0:nlongs+1,1:nlevs),   &
                  ModelVar % u(1:nlongs,1:nlevs),        &
                  ModelVar % v(1:nlongs,1:nlevs))
  CALL Boundaries (ModelVar, set_u=.TRUE., set_v=.TRUE.)


  ! ----------------------------------------------
  IF ((option_gb == 1) .OR. (option_gb == 2)) THEN
    ! 2. Compute the balanced r (r_b) from psi (not known up to a constant on each level)
    IF (option_gb ==1) THEN
      ! Analytical balance
      CALL LinearBal_r (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                        r_b(1:nlongs, 1:nlevs))
    ELSE
      ! Statistical balance
    END IF

    IF (option_reg == 1) THEN
      ! 3. Perform vertical regression
      DO x = 1, nlongs
        r_b(x,1:nlevs) = MATMUL(Regression(1:nlevs, 1:nlevs), r_b(x,1:nlevs))
      END DO
    END IF

    ! 4. Compute the total r
    ModelVar % r(1:nlongs,1:nlevs) = r_b(1:nlongs,1:nlevs) + ControlVar % v3(1:nlongs,1:nlevs)

  ELSE
    ! No geostrophic balance relations used
    ModelVar % r(1:nlongs,1:nlevs) = ControlVar % v3(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries (ModelVar, set_r=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  IF ((option_hb == 1) .OR. (option_hb == 2)) THEN
    IF (option_hb == 1) THEN
      ! 5. Compute the balanced b - see Eq (19) of model paper
      CALL HydroBal_b (ModelVar % r(1:nlongs,0:nlevs+1),   &
                       b_b(1:nlongs,1:nlevs),              &
                       dims)
    ELSE
      ! Statistical balance
    END IF

    ! 6. Compute the total b
    ModelVar % b(1:nlongs,1:nlevs) = b_b(1:nlongs,1:nlevs) + ControlVar % v4(1:nlongs,1:nlevs)
  ELSE
    ! No hydrostatic balance relations used
    ModelVar % b(1:nlongs,1:nlevs) = ControlVar % v4(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries (ModelVar, set_b=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  IF (option_ab == 1) THEN
    ! 7. Compute the balanced w
    CALL Anbalw (LS % rho(0:nlongs+1,0:nlevs+1),       &
                 ModelVar % u(0:nlongs+1,0:nlevs+1),   &
                 w_b(1:nlongs,1:nlevs),                &
                 dims)

    ! 8. Compute the total w
    ModelVar % w(1:nlongs,1:nlevs) = w_b(1:nlongs,1:nlevs) + ControlVar % v5(1:nlongs,1:nlevs)
  ELSE
    ! No anelastic balance relations used
    ModelVar % w(1:nlongs,1:nlevs) = ControlVar % v5(1:nlongs,1:nlevs)
  END IF
  CALL Boundaries (ModelVar, set_w=.TRUE.)
  ! ----------------------------------------------


  ! ----------------------------------------------
  ! 9. Compute the tracer
  ModelVar % tracer(1:nlongs,1:nlevs) = ControlVar % v6(1:nlongs,1:nlevs)
  CALL Boundaries (ModelVar, set_tracer=.TRUE.)
  ! ----------------------------------------------


ELSE
  ! Control variable transform as REP's thesis
  ! ----------------------------------------------

END IF


END SUBROUTINE U_p
