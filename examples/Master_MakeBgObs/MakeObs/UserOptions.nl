&UserOptions
! Generate observations
! ------------------------------
  Generate_mode   = 2
  datadir_ObsSpec = '../ObsSpecification'
  ObsSpec_file    = 'ObsSpec.dat'
  datadir_Obs     = '.'
  Obs_file        = 'Obs.dat'
  datadirABC_in   = '/home/data'
  init_ABC_file   = 'Truth.nc'
  output_ABC_file = 'Truth.nc'
  dt_da           = 600.0
  t0              = 0
  runlength       = 3600.0
  A               = 0.02
  B               = 0.01
  C               = 1.0E4
/
