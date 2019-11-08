! Define namelist
! ---------------
&UserOptions

! Master_PrepareABC_InitState
! ------------------------------
  Init_ABC_opt             = 1
  datadirUM                = '/home/umdata'
  init_um_file             = 'UM_InitConds.nc'
  datadirABC_out           = '.'
  init_ABC_file            = 'ABC_InitialConds.nc'
  latitude                 = 144
  Regular_vert_grid        = .TRUE.
  Adv_tracer               = .TRUE.
  gravity_wave_switch      = .FALSE.
  f                        = 1.0E-4
  A                        = 0.02
  B                        = 0.005
  C                        = 1.0E5
  BoundSpread              = 100.
/
