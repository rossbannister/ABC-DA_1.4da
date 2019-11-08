! Define namelist
! ---------------
&UserOptions

! Reading and processing UM data
! ------------------------------
  datadirImpliedCov        = '.',
  datadirCVT               = '../Master_Calibration',
  CVT_file                 = 'CVT.nc',
  datadirABCfcs            = '../Master_Calibration/Master_Calibration_stage1',
  LS_file                  = 'FC_Ens001_Item001.nc',
  ImplCov_npoints          = 1
  longindex(1)             = 60
  levindex(1)              = 20
/
