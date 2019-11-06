SUBROUTINE U_p_adj (LS, ControlVar, ModelVar, order,             &
                    option_gb, option_hb, option_ab, option_reg, &
                    Regression, dims)

! Code to perform the parameter cvt: ControlVar = U_p* ModelVar

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  dims_type,              &
  nlongs,                 &
  nlevs

IMPLICIT NONE

INCLUDE "Boundaries_adj.interface"


TYPE(ABC_type),  INTENT(IN)    :: LS
TYPE(CV_type),   INTENT(INOUT) :: ControlVar
TYPE(ABC_type),  INTENT(IN)    :: ModelVar
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
TYPE(ABC_type)                 :: Intermediate


! *******************************************************************
! *** The same numbering system in the comments is used as in U_p ***
! *******************************************************************

  Intermediate % u(0:nlongs+1,0:nlevs+1)      = ModelVar % u(0:nlongs+1,0:nlevs+1)
  Intermediate % v(0:nlongs+1,0:nlevs+1)      = ModelVar % v(0:nlongs+1,0:nlevs+1)
  Intermediate % w(0:nlongs+1,0:nlevs+1)      = ModelVar % w(0:nlongs+1,0:nlevs+1)
  Intermediate % r(0:nlongs+1,0:nlevs+1)      = ModelVar % r(0:nlongs+1,0:nlevs+1)
  Intermediate % b(0:nlongs+1,0:nlevs+1)      = ModelVar % b(0:nlongs+1,0:nlevs+1)
  Intermediate % tracer(0:nlongs+1,0:nlevs+1) = ModelVar % tracer(0:nlongs+1,0:nlevs+1)


IF ((order == 1).OR.(order == 2)) THEN
  ! Traditional kind of control variable transform
  ! ----------------------------------------------

  ! Control variables are psi, chi, (unbalanced) r, (unbalanced) b, (unbalanced) w, tracer


  ! ----------------------------------------------
  ! 9. Compute the tracer
  CALL Boundaries_adj (Intermediate, set_tracer=.TRUE.)
  ControlVar % v6(1:nlongs,1:nlevs) = Intermediate % tracer(1:nlongs,1:nlevs)
  ! ----------------------------------------------


  ! ----------------------------------------------
  CALL Boundaries_adj (Intermediate, set_w=.TRUE.)
  IF (option_ab == 1) THEN
    ! 8. Compute the total w
    w_b(1:nlongs,1:nlevs)             = Intermediate % w(1:nlongs,1:nlevs)
    ControlVar % v5(1:nlongs,1:nlevs) = Intermediate % w(1:nlongs,1:nlevs)

    ! 7. Compute the balanced w
    CALL Anbalw_adj (LS % rho(0:nlongs+1,0:nlevs+1),         &
                     Intermediate % u(0:nlongs+1,0:nlevs+1), &
                     w_b(1:nlongs,1:nlevs),                  &
                     dims)
  ELSE
    ! No anelastic balance relations used
    ControlVar % v5(1:nlongs,1:nlevs) = Intermediate % w(1:nlongs,1:nlevs)
  END IF
  ! ----------------------------------------------


  ! ----------------------------------------------
  CALL Boundaries_adj (Intermediate, set_b=.TRUE.)
  IF ((option_hb == 1) .OR. (option_hb == 2)) THEN
    ! 6. Compute the total b
    b_b(1:nlongs,1:nlevs)             = Intermediate % b(1:nlongs,1:nlevs)
    ControlVar % v4(1:nlongs,1:nlevs) = Intermediate % b(1:nlongs,1:nlevs)

    IF (option_hb == 1) THEN
      ! 5. Compute the balanced b - see Eq (19) of model paper
      CALL HydroBal_b_adj (Intermediate % r(1:nlongs,0:nlevs+1),   &
                           b_b(1:nlongs,1:nlevs),                  &
                           dims)
    ELSE
      ! Statistical balance
    END IF

  ELSE
    ! No hydrostatic balance relations used
    ControlVar % v4(1:nlongs,1:nlevs) = Intermediate % b(1:nlongs,1:nlevs)
  END IF
  ! ----------------------------------------------


  ! ----------------------------------------------
  CALL Boundaries_adj (Intermediate, set_r=.TRUE.)
  IF ((option_gb == 1) .OR. (option_gb == 2)) THEN
    ! 4. Compute the total r
    r_b(1:nlongs,1:nlevs)             = Intermediate % r(1:nlongs,1:nlevs)
    ControlVar % v3(1:nlongs,1:nlevs) = Intermediate % r(1:nlongs,1:nlevs)

    IF (option_reg == 1) THEN
      ! 3. Perform vertical regression
      DO x = 1, nlongs
        r_b(x,1:nlevs) = MATMUL(TRANSPOSE(Regression(1:nlevs, 1:nlevs)), r_b(x,1:nlevs))
      END DO
    END IF

    ! 2. Compute the balanced r (r_b) from psi (not known up to a constant on each level)
    IF (option_gb ==1) THEN
      ! Analytical balance
      CALL LinearBal_r_adj (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                            r_b(1:nlongs, 1:nlevs))
    ELSE
      ! Statistical balance
    END IF

  ELSE
    ! No geostrophic balance relations used
    ControlVar % v3(1:nlongs,1:nlevs) = Intermediate % r(1:nlongs,1:nlevs)
  END IF
  ! ----------------------------------------------



  ! ----------------------------------------------
  ! 1. Compute u and v from psi and chi
  ! ----------------------------------------------
  CALL Boundaries_adj (Intermediate, set_u=.TRUE., set_v=.TRUE.)
  CALL Helmholtz_adj (ControlVar % v1(0:nlongs+1,1:nlevs),   &
                      ControlVar % v2(0:nlongs+1,1:nlevs),   &
                      Intermediate % u(1:nlongs,1:nlevs),    &
                      Intermediate % v(1:nlongs,1:nlevs))


ELSE
  ! Control variable transform as REP's thesis
  ! ----------------------------------------------

END IF




END SUBROUTINE U_p_adj
