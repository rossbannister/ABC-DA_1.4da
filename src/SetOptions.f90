SUBROUTINE SetOptions

! Read namelists and set derived parameters

USE DefConsTypes, ONLY :     &
  nlongs, nlevs, pi,         &
  zero,                      &
  sr_nlongs,                 &
  half_sr_nlongs,            &
  UserOptions,               &
  datadirUM,                 &
  Init_ABC_opt,              &
  init_um_file,              &
  latitude,                  &
  Regular_vert_grid,         &
  gravity_wave_switch,       &
  datadirABC_in,             &
  datadirABC_out,            &
  init_ABC_file,             &
  output_ABC_file,           &
  diagnostics_file,          &
  dt, dx, H,                 &
  runlength,                 &
  ndumps, ntimesteps,        &
  f, A, B, C,                &
  BoundSpread,               &
  Adv_tracer,                &
  Lengthscale_diagnostics,   &
  convection_switch,         &
  press_source_x,            &
  press_source_z,            &
  x_scale, z_scale,          &
  press_amp,                 &
  Lx, deltat, dz,            &
  recipdx2, recipdx,         &
  recip2dx,                  &
  third, half,               &
  recippi,                   &
  fourpi2,                   &
  alpha_f, alpha_N,          &
  beta_f, beta_N,            &
  recip_alpha_f,             &
  recip_alpha_N,             &
  bdiva_f, bdiva_N,          &
  CVT_file,                  &
  datadirCVT,                &
  datadirLinearAnal,         &
  datadirTestDA,             &
  RunAdjTests_CVT,           &
  RunAdjTests_obs,           &
  RunInvTests,               &
  LS_file,                   &
  Pert_file,                 &
  CalibRunStage,             &
  NEns,                      &
  EnsDirs,                   &
  NEnsMems,                  &
  NNMC,                      &
  NMCDirs,                   &
  Nlats,                     &
  latindex,                  &
  datadirABCfcs,             &
  datadirABCperts,           &
  datadirRegression,         &
  datadirConParams,          &
  CVT,                       &
  VertSmoothPoints,          &
  HorizSmoothPoints,         &
  ForceCor,                  &
  Generate_mode,             &
  ObsSpec,                   &
  datadir_ObsSpec,           &
  ObsSpec_file,              &
  datadir_Bg,                &
  datadir_Obs,               &
  Bg_file,                   &
  Obs_file,                  &
  dt_da,                     &
  t0,                        &
  random_seed,               &
  datadirImpliedCov,         &
  ImplCov_npoints,           &
  longindex,                 &
  levindex,                  &
  datadirRawCov,             &
  Hybrid_opt,                &
  Vartype,                   &
  datadirAnal,               &
  anal_file,                 &
  analinc_file,              &
  N_outerloops,              &
  N_innerloops_max,          &
  mu, minus_mu,              &
  crit_inner




IMPLICIT NONE

INTEGER                   :: unitno, ErrStatus
INTEGER                   :: dirno, lat, point


! Set the default values for calibration
CALL Initialise_CVT (CVT)
CVT % CVT_order         = 1
CVT % CVT_param_opt_gb  = 1
CVT % CVT_param_opt_hb  = 1
CVT % CVT_param_opt_ab  = 1
CVT % CVT_param_opt_reg = 1
CVT % CVT_vert_opt_sym  = 1
CVT % CVT_stddev_opt    = 2

! Set the default options for the generation of obs and bg
! --------------------------------------------------------
CALL Initialise_ObsSpec (ObsSpec)


! Read the namelist containing user's choice of parameters
! --------------------------------------------------------
ErrStatus = 0
unitno    = 12
OPEN (UNIT=unitno, FILE='UserOptions.nl', IOSTAT=ErrStatus)
IF (ErrStatus == 0) THEN
  READ (unitno, NML=UserOptions)
  CLOSE (unitno)
ELSE
  PRINT*, 'Expecting a namelist file called UserOptions.nl'
  STOP
END IF

! Set the derived parameters
! --------------------------
Lx             = dx * REAL(nlongs)
deltat         = dt / 2.
dz             = H / nlevs
ntimesteps     = runlength / dt
recipdx2       = 1. / (dx * dx)
recipdx        = 1. / dx
recip2dx       = 1. / ( 2. * dx )
third          = 1. / 3.
half           = 1. / 2.
recippi        = 1. / pi
alpha_f        = 1. + deltat*deltat * f*f / 4.
alpha_N        = 1. + deltat*deltat * A*A / 4.
beta_f         = 1. - deltat*deltat * f*f / 4.
beta_N         = 1. - deltat*deltat * A*A / 4.
recip_alpha_f  = 1. / alpha_f
recip_alpha_N  = 1. / alpha_N
bdiva_f        = beta_f / alpha_f
bdiva_N        = beta_N / alpha_N
sr_nlongs      = SQRT(REAL(nlongs))
half_sr_nlongs = sr_nlongs / 2.
fourpi2        = 4. * pi * pi
minus_mu       = -1. * mu

PRINT*
PRINT*, 'pi                       = ', pi
PRINT*, 'zero                     = ', zero
PRINT*
PRINT*, 'User options'
PRINT*, '============================================================'
PRINT*, 'datadirUM                = ', TRIM(datadirUM)
PRINT*, 'init_um_file             = ', TRIM(init_um_file)
PRINT*, 'latitude                 = ', latitude
PRINT*, 'Regular_vert_grid        = ', Regular_vert_grid
PRINT*, 'gravity_wave_switch      = ', gravity_wave_switch
PRINT*, 'Init_ABC_opt             = ', Init_ABC_opt
PRINT*, 'f                        = ', f
PRINT*, 'A                        = ', A
PRINT*, 'B                        = ', B
PRINT*, 'C                        = ', C
PRINT*, 'BoundSpread              = ', BoundSpread
PRINT*, 'datadirABC_in            = ', TRIM(datadirABC_in)
PRINT*, 'datadirABC_out           = ', TRIM(datadirABC_out)
PRINT*, 'init_ABC_file            = ', TRIM(init_ABC_file)
PRINT*, 'output_ABC_file          = ', TRIM(output_ABC_file)
PRINT*, 'diagnostics_file         = ', TRIM(diagnostics_file)
PRINT*, 'dt                       = ', dt
PRINT*, 'dx                       = ', dx
PRINT*, 'H                        = ', H
PRINT*, 'runlength                = ', runlength
PRINT*, 'ndumps                   = ', ndumps
PRINT*, 'convection_switch        = ', convection_switch
PRINT*, 'press_source_x           = ', press_source_x
PRINT*, 'press_source_z           = ', press_source_z
PRINT*, 'x_scale                  = ', x_scale
PRINT*, 'z_scale                  = ', z_scale
PRINT*, 'press_amp                = ', press_amp
PRINT*, 'Adv_tracer               = ', Adv_tracer
PRINT*, 'Lengthscale_diagnostics  = ', Lengthscale_diagnostics
PRINT*, '============================================================'
PRINT*, 'datadirLinearAnal        = ', TRIM(datadirLinearAnal)
PRINT*, '============================================================'
PRINT*, 'datadirTestDA            = ', TRIM(datadirTestDA)
PRINT*, 'RunAdjTests_CVT          = ', RunAdjTests_CVT
PRINT*, 'RunAdjTests_obs          = ', RunAdjTests_obs
PRINT*, 'RunInvTests              = ', RunInvTests
PRINT*, 'LS_file                  = ', TRIM(LS_file)
PRINT*, 'Pert_file                = ', TRIM(Pert_file)
PRINT*, '============================================================'
PRINT*, 'datadirImpliedCov        = ', TRIM(datadirImpliedCov)
PRINT*, 'datadirRawCov            = ', TRIM(datadirRawCov)
PRINT*, 'ImplCov_npoints          = ', ImplCov_npoints
DO point = 1, ImplCov_npoints
  PRINT *, 'Long index', point, ' = ', longindex(point)
  PRINT *, 'Lev index', point, '  = ', levindex(point)
END DO
PRINT*, '============================================================'
PRINT*, 'CalibRunStage            = ', CalibRunStage
PRINT*, 'datadirCVT               = ', TRIM(datadirCVT)
PRINT*, 'CVT_file                 = ', TRIM(CVT_file)
PRINT*, 'NEns                     = ', NEns
PRINT*, 't0                       = ', t0
DO dirno = 1, NEns
  PRINT*, 'EnsDirs', dirno, '= ', TRIM(EnsDirs(dirno))
END DO
PRINT*, 'NEnsMems                 = ', NEnsMems
PRINT*, 'NNMC                     = ', NNMC
DO dirno = 1, NNMC
  PRINT*, 'NCMDirs', dirno, '= ', TRIM(NMCDirs(dirno))
END DO
PRINT*, 'Nlats                    = ', Nlats
DO lat = 1, Nlats
  PRINT*, 'Lat index ', lat, '= ', latindex(lat)
END DO
PRINT*, 'datadirABCfcs            = ', TRIM(datadirABCfcs)
PRINT*, 'datadirABCperts          = ', TRIM(datadirABCperts)
PRINT*, 'datadirRegression        = ', TRIM(datadirRegression)
PRINT*, 'datadirConParams         = ', TRIM(datadirConParams)
PRINT*, 'CVT % CVT_order          = ', CVT % CVT_order
PRINT*, 'CVT % CVT_param_opt_gb   = ', CVT % CVT_param_opt_gb
PRINT*, 'CVT % CVT_param_opt_hb   = ', CVT % CVT_param_opt_hb
PRINT*, 'CVT % CVT_param_opt_ab   = ', CVT % CVT_param_opt_ab
PRINT*, 'CVT % CVT_param_opt_reg  = ', CVT % CVT_param_opt_reg
PRINT*, 'CVT % CVT_vert_opt_sym   = ', CVT % CVT_vert_opt_sym
PRINT*, 'CVT % CVT_stddev_opt     = ', CVT % CVT_stddev_opt
PRINT*, 'VertSmoothPoints         = ', VertSmoothPoints
PRINT*, 'HorizSmoothPoints        = ', HorizSmoothPoints
PRINT*, 'ForceCor                 = ', ForceCor
PRINT*, 'Generate_mode            = ', Generate_mode
PRINT*, 'datadir_ObsSpec          = ', TRIM(datadir_ObsSpec)
PRINT*, 'ObsSpec_file             = ', TRIM(ObsSpec_file)
PRINT*, 'ObsSpec % year0          = ', ObsSpec % year0
PRINT*, 'ObsSpec % month0         = ', ObsSpec % month0
PRINT*, 'ObsSpec % day0           = ', ObsSpec % day0
PRINT*, 'ObsSpec % hour0          = ', ObsSpec % hour0
PRINT*, 'ObsSpec % min0           = ', ObsSpec % min0
PRINT*, 'ObsSpec % sec0           = ', ObsSpec % sec0
PRINT*, 'ObsSpec % NumBatches     = ', ObsSpec % NumBatches
PRINT*, 'ObsSpec % batch()        = ', ObsSpec % batch(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % seconds()      = ', ObsSpec % seconds(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % ob_of_what()   = ', ObsSpec % ob_of_what(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % NumObs_long()  = ', ObsSpec % NumObs_long(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % NumObs_height()= ', ObsSpec % NumObs_height(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % long_min()     = ', ObsSpec % long_min(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % long_max()     = ', ObsSpec % long_max(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % height_min()   = ', ObsSpec % height_min(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % height_max()   = ', ObsSpec % height_max(1:ObsSpec % NumBatches)
PRINT*, 'ObsSpec % stddev()       = ', ObsSpec % stddev(1:ObsSpec % NumBatches)
PRINT*, 'datadir_Bg               = ', TRIM(datadir_Bg)
PRINT*, 'Bg_file                  = ', TRIM(Bg_file)
PRINT*, 'datadir_Obs              = ', TRIM(datadir_Obs)
PRINT*, 'Obs_file                 = ', TRIM(Obs_file)
PRINT*, 'dt_da                    = ', dt_da
PRINT*, '============================================================'
PRINT*, 'datadirAnal              = ', TRIM(datadirAnal)
PRINT*, 'anal_file                = ', TRIM(anal_file)
PRINT*, 'analinc_file             = ', TRIM(analinc_file)
PRINT*, 'Hybrid_opt               = ', Hybrid_opt
PRINT*, 'Vartype                  = ', Vartype
PRINT*, 'N_outerloops             = ', N_outerloops
PRINT*, 'N_innerloops_max         = ', N_innerloops_max
PRINT*, 'crit_inner               = ', crit_inner
PRINT*, 'mu                       = ', mu
PRINT*, 'minus_mu                 = ', minus_mu
PRINT*, 'Derived variables'
PRINT*, '============================================================'
PRINT*, 'Lx                       = ', Lx
PRINT*, 'deltat                   = ', deltat
PRINT*, 'dz                       = ', dz
PRINT*, 'ntimesteps               = ', ntimesteps
PRINT*, 'recipdx2                 = ', recipdx2
PRINT*, 'recipdx                  = ', recipdx
PRINT*, 'third                    = ', third
PRINT*, 'recippi                  = ', recippi
PRINT*, 'alpha_f                  = ', alpha_f
PRINT*, 'beta_f                   = ', beta_f
PRINT*, 'beta_N                   = ', beta_N
PRINT*, 'recip_alpha_f            = ', recip_alpha_f
PRINT*, 'recip_alpha_N            = ', recip_alpha_N
PRINT*, 'bdiva_f                  = ', bdiva_f
PRINT*, 'bdiva_N                  = ', bdiva_N
PRINT*, 'sr_nlongs                = ', sr_nlongs
PRINT*, 'half_sr_nlongs           = ', half_sr_nlongs
PRINT*, 'random_seed              = ', random_seed
PRINT*, '============================================================'


! Seed the random number generator
PRINT*, 'Seeding the random number generator'
CALL SRAND(random_seed)


END SUBROUTINE SetOptions
