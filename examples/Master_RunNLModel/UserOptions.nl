! Define namelist
! ---------------
&UserOptions

! Running ABC model
! ------------------------------
  datadirABC_in            = '../Master_PrepareABC_InitState'
  init_ABC_file            = 'ABC_InitialConds.nc'
  datadirABC_out           = '.'
  output_ABC_file          = 'ABC_ModelRun.nc'
  diagnostics_file         = 'ABC_Diagnostics.dat'
  runlength                = 10800.0
  ndumps                   = 6
  dt                       = 1.0
  Lengthscale_diagnostics  = .TRUE.
  A                        = 0.02
  B                        = 0.01
  C                        = 1.0E4
  Adv_tracer               = .TRUE.
/
