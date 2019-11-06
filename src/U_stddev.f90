SUBROUTINE U_stddev (ControlVar,                                      &
                     CVT)

! Code to multiply by each control parameter's standard deviation field

USE DefConsTypes, ONLY :  &
  ZREAL8,                 &
  ABC_type,               &
  CV_type,                &
  CVT_type,               &
  nlongs,                 &
  nlevs


IMPLICIT NONE

TYPE(CV_type),  INTENT(INOUT) :: ControlVar
TYPE(CVT_type), INTENT(IN)    :: CVT


ControlVar % v1(1:nlongs, 1:nlevs) = ControlVar % v1(1:nlongs, 1:nlevs) * &
                                     CVT % sigma1(1:nlongs, 1:nlevs)
ControlVar % v2(1:nlongs, 1:nlevs) = ControlVar % v2(1:nlongs, 1:nlevs) * &
                                     CVT % sigma2(1:nlongs, 1:nlevs)
ControlVar % v3(1:nlongs, 1:nlevs) = ControlVar % v3(1:nlongs, 1:nlevs) * &
                                     CVT % sigma3(1:nlongs, 1:nlevs)
ControlVar % v4(1:nlongs, 1:nlevs) = ControlVar % v4(1:nlongs, 1:nlevs) * &
                                     CVT % sigma4(1:nlongs, 1:nlevs)
ControlVar % v5(1:nlongs, 1:nlevs) = ControlVar % v5(1:nlongs, 1:nlevs) * &
                                     CVT % sigma5(1:nlongs, 1:nlevs)
ControlVar % v6(1:nlongs, 1:nlevs) = ControlVar % v6(1:nlongs, 1:nlevs) * &
                                     CVT % sigma6(1:nlongs, 1:nlevs)

END SUBROUTINE U_stddev
