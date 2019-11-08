! Define namelist
! ---------------
&UserOptions

! Reading and processing UM data
! ------------------------------
  datadirTestDA            = '.',
  datadirABCfcs            = '../../Master_Calibration/Master_Calibration_stage1',
  LS_file                  = 'FC_Ens001_Item001.nc',
  datadirABCperts          = '../../Master_Calibration/Master_Calibration_stage2
',
  Pert_file                = 'PertABC_Ens001_Item001.nc',
  datadir_Obs              = '../../Master_MakeBgObs/MakeObs',
  Obs_file                 = 'Obs.dat',
  RunAdjTests_CVT          = .FALSE.,
  RunAdjTests_obs          = .TRUE.,
  RunInvTests              = .FALSE.,
  diagnostics_file         = 'Diagnostics.dat'
/
