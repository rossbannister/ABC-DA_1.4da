! Define namelist
! ---------------
&UserOptions

! Run the data assimilation
! -------------------------
  Vartype          = 3,                  !3=3DVar, 35=3D-FGAT, 4=4DVar
  Hybrid_opt       = 1,                  !Type of hybrid (or if pure Var)
  datadir_Bg       = '../Master_MakeBgObs/MakeBg',
  Bg_file          = 'Bg.nc',
  datadirCVT       = '../Master_Calibration',
  CVT_file         = 'CVT.nc',
  datadir_Obs      = '../Master_MakeBgObs/MakeObs',
  Obs_file         = 'Obs.dat',
  t0               = 0,                  !Time of start of this DA cycle
  N_outerloops     = 2,
  N_innerloops_max = 20,
  crit_inner       = 0.01,
  datadirAnal      = '.',
  anal_file        = 'Anal.nc',
  analinc_file     = 'AnalInc.nc',
  diagnostics_file = 'diagnostics.dat'
/
