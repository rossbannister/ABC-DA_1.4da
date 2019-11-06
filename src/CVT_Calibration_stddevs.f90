SUBROUTINE CVT_Calibration_stddevs (CVT, pert, div_cons)

! If div_cons = 0.0, use the perturbation to contribute to the sigma values

! If div_cons > 0.0, square-root and divide sigmas by div_cons
!                    (to compute std dev)

USE DefConsTypes, ONLY :     &
    ZREAL8,                  &
    CVT_type,                &
    CV_type,                 &
    nlongs, nlevs,           &
    small, ConditionFudge

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT
TYPE(CV_type),  INTENT(IN)      :: pert
REAL(ZREAL8),   INTENT(IN)      :: div_cons

REAL(ZREAL8)                    :: mul_cons, mean, nlongs_r, nlevs_r, increment
INTEGER                         :: x, z

IF (DABS(div_cons) < small) THEN
  ! Contribute to the sigma values for this perturbation

  CVT % sigma1(1:nlongs, 1:nlevs) = CVT % sigma1(1:nlongs, 1:nlevs) + &
                                    pert % v1(1:nlongs, 1:nlevs) *    &
                                    pert % v1(1:nlongs, 1:nlevs)
  CVT % sigma2(1:nlongs, 1:nlevs) = CVT % sigma2(1:nlongs, 1:nlevs) + &
                                    pert % v2(1:nlongs, 1:nlevs) *    &
                                    pert % v2(1:nlongs, 1:nlevs)
  CVT % sigma3(1:nlongs, 1:nlevs) = CVT % sigma3(1:nlongs, 1:nlevs) + &
                                    pert % v3(1:nlongs, 1:nlevs) *    &
                                    pert % v3(1:nlongs, 1:nlevs)
  CVT % sigma4(1:nlongs, 1:nlevs) = CVT % sigma4(1:nlongs, 1:nlevs) + &
                                    pert % v4(1:nlongs, 1:nlevs) *    &
                                    pert % v4(1:nlongs, 1:nlevs)
  CVT % sigma5(1:nlongs, 1:nlevs) = CVT % sigma5(1:nlongs, 1:nlevs) + &
                                    pert % v5(1:nlongs, 1:nlevs) *    &
                                    pert % v5(1:nlongs, 1:nlevs)
  CVT % sigma6(1:nlongs, 1:nlevs) = CVT % sigma6(1:nlongs, 1:nlevs) + &
                                    pert % v6(1:nlongs, 1:nlevs) *    &
                                    pert % v6(1:nlongs, 1:nlevs)
ELSE


  ! Normalize

  mul_cons = 1.0 / div_cons
  CVT % sigma1(1:nlongs, 1:nlevs) = CVT % sigma1(1:nlongs, 1:nlevs) * mul_cons
  CVT % sigma2(1:nlongs, 1:nlevs) = CVT % sigma2(1:nlongs, 1:nlevs) * mul_cons
  CVT % sigma3(1:nlongs, 1:nlevs) = CVT % sigma3(1:nlongs, 1:nlevs) * mul_cons
  CVT % sigma4(1:nlongs, 1:nlevs) = CVT % sigma4(1:nlongs, 1:nlevs) * mul_cons
  CVT % sigma5(1:nlongs, 1:nlevs) = CVT % sigma5(1:nlongs, 1:nlevs) * mul_cons
  CVT % sigma6(1:nlongs, 1:nlevs) = CVT % sigma6(1:nlongs, 1:nlevs) * mul_cons

  ! Condition for zero values
  increment = SUM(CVT % sigma1(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma1(1:nlongs, 1:nlevs) = CVT % sigma1(1:nlongs, 1:nlevs) + increment

  increment = SUM(CVT % sigma2(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma2(1:nlongs, 1:nlevs) = CVT % sigma2(1:nlongs, 1:nlevs) + increment

  increment = SUM(CVT % sigma3(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma3(1:nlongs, 1:nlevs) = CVT % sigma3(1:nlongs, 1:nlevs) + increment

  increment = SUM(CVT % sigma4(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma4(1:nlongs, 1:nlevs) = CVT % sigma4(1:nlongs, 1:nlevs) + increment

  increment = SUM(CVT % sigma5(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma5(1:nlongs, 1:nlevs) = CVT % sigma5(1:nlongs, 1:nlevs) + increment

  increment = SUM(CVT % sigma6(1:nlongs, 1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % sigma6(1:nlongs, 1:nlevs) = CVT % sigma6(1:nlongs, 1:nlevs) + increment

  ! Square-root
  CVT % sigma1(1:nlongs, 1:nlevs) = SQRT(CVT % sigma1(1:nlongs, 1:nlevs))
  CVT % sigma2(1:nlongs, 1:nlevs) = SQRT(CVT % sigma2(1:nlongs, 1:nlevs))
  CVT % sigma3(1:nlongs, 1:nlevs) = SQRT(CVT % sigma3(1:nlongs, 1:nlevs))
  CVT % sigma4(1:nlongs, 1:nlevs) = SQRT(CVT % sigma4(1:nlongs, 1:nlevs))
  CVT % sigma5(1:nlongs, 1:nlevs) = SQRT(CVT % sigma5(1:nlongs, 1:nlevs))
  CVT % sigma6(1:nlongs, 1:nlevs) = SQRT(CVT % sigma6(1:nlongs, 1:nlevs))



  ! Do requested averaging according to option
  SELECT CASE (CVT % CVT_stddev_opt)

  CASE (2)  ! Standard deviations are level dependent only
    nlongs_r = REAL(nlongs)
    mul_cons = 1.0 / nlongs_r
    DO z = 1, nlevs
      mean = SUM(CVT % sigma1(1:nlongs,z)) * mul_cons
      CVT % sigma1(1:nlongs,z) = mean
      mean = SUM(CVT % sigma2(1:nlongs,z)) * mul_cons
      CVT % sigma2(1:nlongs,z) = mean
      mean = SUM(CVT % sigma3(1:nlongs,z)) * mul_cons
      CVT % sigma3(1:nlongs,z) = mean
      mean = SUM(CVT % sigma4(1:nlongs,z)) * mul_cons
      CVT % sigma4(1:nlongs,z) = mean
      mean = SUM(CVT % sigma5(1:nlongs,z)) * mul_cons
      CVT % sigma5(1:nlongs,z) = mean
      mean = SUM(CVT % sigma6(1:nlongs,z)) * mul_cons
      CVT % sigma6(1:nlongs,z) = mean
    END DO

  CASE (3)  ! Standard deviations are constant
    nlongs_r = REAL(nlongs)
    nlevs_r  = REAL(nlevs)
    mul_cons = 1.0 / (nlongs_r * nlevs_r)
    mean = SUM(CVT % sigma1(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma1(1:nlongs,1:nlevs) = mean
    mean = SUM(CVT % sigma2(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma2(1:nlongs,1:nlevs) = mean
    mean = SUM(CVT % sigma3(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma3(1:nlongs,1:nlevs) = mean
    mean = SUM(CVT % sigma4(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma4(1:nlongs,1:nlevs) = mean
    mean = SUM(CVT % sigma5(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma5(1:nlongs,1:nlevs) = mean
    mean = SUM(CVT % sigma6(1:nlongs,1:nlevs)) * mul_cons
    CVT % sigma6(1:nlongs,1:nlevs) = mean

  END SELECT


END IF

END SUBROUTINE CVT_Calibration_stddevs
