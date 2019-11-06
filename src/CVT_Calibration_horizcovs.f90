SUBROUTINE CVT_Calibration_horizcovs (CVT, pert, div_cons)

! If div_cons = 0.0, use the perturbation to contribute to the horizontal variance spectra

! If div_cons > 0.0, divide horiz covs by div_cons
!                    (to compute horizontal variance spectra)

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

REAL(ZREAL8)                    :: mul_cons, increment
INTEGER                         :: k
INTEGER                         :: real_index, imag_index


IF (DABS(div_cons) < small) THEN
  ! Contribute to the horizontal variance spectra for this perturbation

  ! Deal with the largest scale
  CVT % HorizEV1(1,1:nlevs) = CVT % HorizEV1(1,1:nlevs) + &
                              pert % v1(1,1:nlevs) * pert % v1(1,1:nlevs)
  CVT % HorizEV2(1,1:nlevs) = CVT % HorizEV2(1,1:nlevs) + &
                              pert % v2(1,1:nlevs) * pert % v2(1,1:nlevs)
  CVT % HorizEV3(1,1:nlevs) = CVT % HorizEV3(1,1:nlevs) + &
                              pert % v3(1,1:nlevs) * pert % v3(1,1:nlevs)
  CVT % HorizEV4(1,1:nlevs) = CVT % HorizEV4(1,1:nlevs) + &
                              pert % v4(1,1:nlevs) * pert % v4(1,1:nlevs)
  CVT % HorizEV5(1,1:nlevs) = CVT % HorizEV5(1,1:nlevs) + &
                              pert % v5(1,1:nlevs) * pert % v5(1,1:nlevs)
  CVT % HorizEV6(1,1:nlevs) = CVT % HorizEV6(1,1:nlevs) + &
                              pert % v6(1,1:nlevs) * pert % v6(1,1:nlevs)

  ! Deal with the intermediate scales
  DO k = 2, nlongs/2
    real_index = 2*k-2
    imag_index = 2*k-1
    CVT % HorizEV1(k,1:nlevs) = CVT % HorizEV1(k,1:nlevs) +                                     &
                                pert % v1(real_index,1:nlevs) * pert % v1(real_index,1:nlevs) + &
                                pert % v1(imag_index,1:nlevs) * pert % v1(imag_index,1:nlevs)
    CVT % HorizEV2(k,1:nlevs) = CVT % HorizEV2(k,1:nlevs) +                                     &
                                pert % v2(real_index,1:nlevs) * pert % v2(real_index,1:nlevs) + &
                                pert % v2(imag_index,1:nlevs) * pert % v2(imag_index,1:nlevs)
    CVT % HorizEV3(k,1:nlevs) = CVT % HorizEV3(k,1:nlevs) +                                     &
                                pert % v3(real_index,1:nlevs) * pert % v3(real_index,1:nlevs) + &
                                pert % v3(imag_index,1:nlevs) * pert % v3(imag_index,1:nlevs)
    CVT % HorizEV4(k,1:nlevs) = CVT % HorizEV4(k,1:nlevs) +                                     &
                                pert % v4(real_index,1:nlevs) * pert % v4(real_index,1:nlevs) + &
                                pert % v4(imag_index,1:nlevs) * pert % v4(imag_index,1:nlevs)
    CVT % HorizEV5(k,1:nlevs) = CVT % HorizEV5(k,1:nlevs) +                                     &
                                pert % v5(real_index,1:nlevs) * pert % v5(real_index,1:nlevs) + &
                                pert % v5(imag_index,1:nlevs) * pert % v5(imag_index,1:nlevs)
    CVT % HorizEV6(k,1:nlevs) = CVT % HorizEV6(k,1:nlevs) +                                     &
                                pert % v6(real_index,1:nlevs) * pert % v6(real_index,1:nlevs) + &
                                pert % v6(imag_index,1:nlevs) * pert % v6(imag_index,1:nlevs)
  END DO

  ! Deal with the smallest scale
  CVT % HorizEV1(nlongs/2+1,1:nlevs) = CVT % HorizEV1(nlongs/2+1,1:nlevs) + &
                                       pert % v1(nlongs,1:nlevs) * pert % v1(nlongs,1:nlevs)
  CVT % HorizEV2(nlongs/2+1,1:nlevs) = CVT % HorizEV2(nlongs/2+1,1:nlevs) + &
                                       pert % v2(nlongs,1:nlevs) * pert % v2(nlongs,1:nlevs)
  CVT % HorizEV3(nlongs/2+1,1:nlevs) = CVT % HorizEV3(nlongs/2+1,1:nlevs) + &
                                       pert % v3(nlongs,1:nlevs) * pert % v3(nlongs,1:nlevs)
  CVT % HorizEV4(nlongs/2+1,1:nlevs) = CVT % HorizEV4(nlongs/2+1,1:nlevs) + &
                                       pert % v4(nlongs,1:nlevs) * pert % v4(nlongs,1:nlevs)
  CVT % HorizEV5(nlongs/2+1,1:nlevs) = CVT % HorizEV5(nlongs/2+1,1:nlevs) + &
                                       pert % v5(nlongs,1:nlevs) * pert % v5(nlongs,1:nlevs)
  CVT % HorizEV6(nlongs/2+1,1:nlevs) = CVT % HorizEV6(nlongs/2+1,1:nlevs) + &
                                       pert % v6(nlongs,1:nlevs) * pert % v6(nlongs,1:nlevs)


ELSE


  ! Normalize
  mul_cons = 1.0 / div_cons
  CVT % HorizEV1(1:nlongs/2+1,1:nlevs) = CVT % HorizEV1(1:nlongs/2+1,1:nlevs) * mul_cons
  CVT % HorizEV2(1:nlongs/2+1,1:nlevs) = CVT % HorizEV2(1:nlongs/2+1,1:nlevs) * mul_cons
  CVT % HorizEV3(1:nlongs/2+1,1:nlevs) = CVT % HorizEV3(1:nlongs/2+1,1:nlevs) * mul_cons
  CVT % HorizEV4(1:nlongs/2+1,1:nlevs) = CVT % HorizEV4(1:nlongs/2+1,1:nlevs) * mul_cons
  CVT % HorizEV5(1:nlongs/2+1,1:nlevs) = CVT % HorizEV5(1:nlongs/2+1,1:nlevs) * mul_cons
  CVT % HorizEV6(1:nlongs/2+1,1:nlevs) = CVT % HorizEV6(1:nlongs/2+1,1:nlevs) * mul_cons

  ! Condition for small values
  increment = SUM(CVT % HorizEV1(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV1(1:nlongs/2+1,1:nlevs) = CVT % HorizEV1(1:nlongs/2+1,1:nlevs) + increment

  increment = SUM(CVT % HorizEV2(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV2(1:nlongs/2+1,1:nlevs) = CVT % HorizEV2(1:nlongs/2+1,1:nlevs) + increment

  increment = SUM(CVT % HorizEV3(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV3(1:nlongs/2+1,1:nlevs) = CVT % HorizEV3(1:nlongs/2+1,1:nlevs) + increment

  increment = SUM(CVT % HorizEV4(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV4(1:nlongs/2+1,1:nlevs) = CVT % HorizEV4(1:nlongs/2+1,1:nlevs) + increment

  increment = SUM(CVT % HorizEV5(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV5(1:nlongs/2+1,1:nlevs) = CVT % HorizEV5(1:nlongs/2+1,1:nlevs) + increment

  increment = SUM(CVT % HorizEV6(1:nlongs/2+1,1:nlevs)) * ConditionFudge
  IF (increment < small) increment = small
  CVT % HorizEV6(1:nlongs/2+1,1:nlevs) = CVT % HorizEV6(1:nlongs/2+1,1:nlevs) + increment

  ! Square root
  CVT % HorizEV1(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV1(1:nlongs/2+1,1:nlevs))
  CVT % HorizEV2(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV2(1:nlongs/2+1,1:nlevs))
  CVT % HorizEV3(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV3(1:nlongs/2+1,1:nlevs))
  CVT % HorizEV4(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV4(1:nlongs/2+1,1:nlevs))
  CVT % HorizEV5(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV5(1:nlongs/2+1,1:nlevs))
  CVT % HorizEV6(1:nlongs/2+1,1:nlevs) = SQRT(CVT % HorizEV6(1:nlongs/2+1,1:nlevs))

END IF

END SUBROUTINE CVT_Calibration_horizcovs
