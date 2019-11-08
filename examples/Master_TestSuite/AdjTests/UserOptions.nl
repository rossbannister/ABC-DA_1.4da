! Define namelist
! ---------------
&UserOptions

! Reading and processing UM data
! ------------------------------
  datadirTestDA            = '.',
  datadirCVT               = '../../Master_Calibration',
  CVT_file                 = 'CVT.nc',
  datadirABCfcs            = '../../Master_Calibration/Master_Calibration_stage1',
  LS_file                  = 'FC_Ens001_Item001.nc',
  datadirABCperts          = '../../Master_Calibration/Master_Calibration_stage2',
  Pert_file                = 'PertABC_Ens001_Item001.nc',
  RunAdjTests              = .TRUE.,
  RunInvTests              = .FALSE.,
  diagnostics_file         = 'Diagnostics.dat'
/
