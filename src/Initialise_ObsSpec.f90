SUBROUTINE Initialise_ObsSpec (ObsSpec)

! Initialise variable that specifies observations


USE DefConsTypes, ONLY :     &
    ObsSpec_type,            &
    maxbatches

IMPLICIT NONE

TYPE(ObsSpec_type), INTENT(INOUT)   :: ObsSpec

  ObsSpec % year0                       = 2000
  ObsSpec % month0                      = 1
  ObsSpec % day0                        = 1
  ObsSpec % hour0                       = 0
  ObsSpec % min0                        = 0
  ObsSpec % sec0                        = 0
  ObsSpec % NumBatches                  = 0
  ObsSpec % batch(1:maxbatches)         = 0
  ObsSpec % seconds(1:maxbatches)       = 0
  ObsSpec % ob_of_what(1:maxbatches)    = 0
  ObsSpec % NumObs_long(1:maxbatches)   = 0
  ObsSpec % NumObs_height(1:maxbatches) = 0
  ObsSpec % long_min(1:maxbatches)      = 0.0
  ObsSpec % long_max(1:maxbatches)      = 0.0
  ObsSpec % height_min(1:maxbatches)    = 0.0
  ObsSpec % height_max(1:maxbatches)    = 0.0
  ObsSpec % stddev(1:maxbatches)        = 0.0

END SUBROUTINE Initialise_ObsSpec
