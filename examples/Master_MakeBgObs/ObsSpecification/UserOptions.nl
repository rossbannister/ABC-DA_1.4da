! Define namelist
! ---------------
&UserOptions

! Specification of observations
! ------------------------------
  Generate_mode            = 1
  ObsSpec%year0            = 2010
  ObsSpec%month0           = 1
  ObsSpec%day0             = 1
  ObsSpec%hour0            = 0
  ObsSpec%min0             = 0
  ObsSpec%sec0             = 0
  ObsSpec%NumBatches       = 8
  datadir_ObsSpec          = '.'
  ObsSpec_file             = 'ObsSpec.dat'
! ------------------------------
  ObsSpec%batch(1)         = 1
  ObsSpec%seconds(1)       = 300
  ObsSpec%ob_of_what(1)    = 1       ! u
  ObsSpec%NumObs_long(1)   = 2
  ObsSpec%NumObs_height(1) = 2
  ObsSpec%long_min(1)      = 10000.0
  ObsSpec%long_max(1)      = 260000.0
  ObsSpec%height_min(1)    = 1000.0
  ObsSpec%height_max(1)    = 12000.0
  ObsSpec%stddev(1)        = 0.5
! ------------------------------
  ObsSpec%batch(2)         = 1
  ObsSpec%seconds(2)       = 600
  ObsSpec%ob_of_what(2)    = 2       ! v
  ObsSpec%NumObs_long(2)   = 2
  ObsSpec%NumObs_height(2) = 2
  ObsSpec%long_min(2)      = 50000.0
  ObsSpec%long_max(2)      = 350000.0
  ObsSpec%height_min(2)    = 1000.0
  ObsSpec%height_max(2)    = 6000.0
  ObsSpec%stddev(2)        = 0.5
! ------------------------------
  ObsSpec%batch(3)         = 1
  ObsSpec%seconds(3)       = 900
  ObsSpec%ob_of_what(3)    = 3       ! w
  ObsSpec%NumObs_long(3)   = 2
  ObsSpec%NumObs_height(3) = 2
  ObsSpec%long_min(3)      = 200000.0
  ObsSpec%long_max(3)      = 400000.0
  ObsSpec%height_min(3)    = 2000.0
  ObsSpec%height_max(3)    = 10000.0
  ObsSpec%stddev(3)        = 0.5
! ------------------------------
  ObsSpec%batch(4)         = 1
  ObsSpec%seconds(4)       = 1200
  ObsSpec%ob_of_what(4)    = 4       ! r
  ObsSpec%NumObs_long(4)   = 2
  ObsSpec%NumObs_height(4) = 2
  ObsSpec%long_min(4)      = 50000.0
  ObsSpec%long_max(4)      = 350000.0
  ObsSpec%height_min(4)    = 1000.0
  ObsSpec%height_max(4)    = 6000.0
  ObsSpec%stddev(4)        = 0.01
! ------------------------------
  ObsSpec%batch(5)         = 1
  ObsSpec%seconds(5)       = 1500
  ObsSpec%ob_of_what(5)    = 5       ! b
  ObsSpec%NumObs_long(5)   = 2
  ObsSpec%NumObs_height(5) = 2
  ObsSpec%long_min(5)      = 1000.0
  ObsSpec%long_max(5)      = 30000.0
  ObsSpec%height_min(5)    = 500.0
  ObsSpec%height_max(5)    = 1000.0
  ObsSpec%stddev(5)        = 0.0001
! ------------------------------
  ObsSpec%batch(6)         = 1
  ObsSpec%seconds(6)       = 1800
  ObsSpec%ob_of_what(6)    = 6       ! tracer
  ObsSpec%NumObs_long(6)   = 2
  ObsSpec%NumObs_height(6) = 2
  ObsSpec%long_min(6)      = 100000.0
  ObsSpec%long_max(6)      = 500000.0
  ObsSpec%height_min(6)    = 20000.0
  ObsSpec%height_max(6)    = 30000.0
  ObsSpec%stddev(6)        = 0.1
! ------------------------------
  ObsSpec%batch(7)         = 1
  ObsSpec%seconds(7)       = 2100
  ObsSpec%ob_of_what(7)    = 7       ! horizontal wind speed
  ObsSpec%NumObs_long(7)   = 2
  ObsSpec%NumObs_height(7) = 2
  ObsSpec%long_min(7)      = 500000.0
  ObsSpec%long_max(7)      = 600000.0
  ObsSpec%height_min(7)    = 2000.0
  ObsSpec%height_max(7)    = 4000.0
  ObsSpec%stddev(7)        = 0.4
! ------------------------------
  ObsSpec%batch(8)         = 1
  ObsSpec%seconds(8)       = 2400
  ObsSpec%ob_of_what(8)    = 8       ! total wind speed
  ObsSpec%NumObs_long(8)   = 2
  ObsSpec%NumObs_height(8) = 2
  ObsSpec%long_min(8)      = 50000.0
  ObsSpec%long_max(8)      = 350000.0
  ObsSpec%height_min(8)    = 1000.0
  ObsSpec%height_max(8)    = 6000.0
  ObsSpec%stddev(8)        = 0.4
/
