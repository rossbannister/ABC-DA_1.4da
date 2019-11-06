!===================================================================================================
MODULE DefConsTypes

!*****************************************************
!*   Definition of all global constants              *
!*   and derived variable types                      *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*   R. Bannister, 1.4da 20-10-2017                  *
!*                                                   *
!*****************************************************

IMPLICIT NONE

!Definition of some new data types
!Define an integer that describes a unified double precision
!-----------------------------------------------------------
INTEGER, PARAMETER        :: ZREAL8         = SELECTED_REAL_KIND(15,307)
INTEGER, PARAMETER        :: wp             = KIND(1.0D0)
INTEGER, PARAMETER        :: Nensmax        = 50
INTEGER, PARAMETER        :: NNMCmax        = 50
INTEGER, PARAMETER        :: Nlatsmax       = 260
INTEGER, PARAMETER        :: Npointsmax     = 24
INTEGER, PARAMETER        :: maxbatches     = 100
REAL(ZREAL8)              :: ConditionFudge = 0.000000001
REAL(ZREAL8)              :: small          = 0.00000000001
REAL(ZREAL8)              :: zero           = 0.0
REAL(ZREAL8)              :: unity          = 1.0
INTEGER                   :: random_seed    = 0

! Model parameters
!-----------------
INTEGER, PARAMETER        :: nlongs = 360            ! total number of longitude points
INTEGER, PARAMETER        :: nlevs  = 60             ! total number of vertical levels
INTEGER                   :: ntimesteps              ! number of timesteps to integrate
REAL(ZREAL8)              :: dt     = 4.             ! full timestep
REAL(ZREAL8)              :: dx     = 1500.          ! dx grid resolution
REAL(ZREAL8)              :: H      = 14862.01       ! UM model height
REAL(ZREAL8)              :: deltat
REAL(ZREAL8)              :: Lx
REAL(ZREAL8)              :: dz

! Mathematical and physical constants
!------------------------------------
REAL(ZREAL8), PARAMETER   :: Rd    = 287.058         ! Gas constant for dry air
REAL(ZREAL8), PARAMETER   :: Re    = 6.371E6         ! Mean radius of the earth
REAL(ZREAL8), PARAMETER   :: Cp    = 1005.7          ! Specific heat at constant pressure
REAL(ZREAL8), PARAMETER   :: Cv    = 719.0           ! Specific heat at constant volume
REAL(ZREAL8), PARAMETER   :: g     = 9.81            ! Acceleration due to gravity
REAL(ZREAL8), PARAMETER   :: p00   = 1.0E5           ! Reference surface pressure 1000 hPa
REAL(ZREAL8), PARAMETER   :: kappa = 0.286           ! Rd/Cp
REAL(ZREAL8), PARAMETER   :: pi    = 3.141592654     ! pi

! Tuneable model parameters
!--------------------------
REAL(ZREAL8)              :: f     = 1.0E-4          ! Coriolis Parameter
REAL(ZREAL8)              :: A     = 0.02            ! A is the buoyancy frequency
REAL(ZREAL8)              :: B     = 0.01            ! B accoustic wave speed modulator
REAL(ZREAL8)              :: C     = 1.0E5           ! Constant relating pressure and density
REAL(ZREAL8)              :: BoundSpread = 50.       ! No of grid points to spread boundary discontinuity info
REAL(ZREAL8)              :: theta_r     = 273.

! Useful constants
!-----------------
REAL(ZREAL8)              :: third, half
REAL(ZREAL8)              :: recippi
REAL(ZREAL8)              :: recipdx
REAL(ZREAL8)              :: recipdx2
REAL(ZREAL8)              :: recip2dx
REAL(ZREAL8)              :: alpha_f
REAL(ZREAL8)              :: alpha_N
REAL(ZREAL8)              :: beta_f
REAL(ZREAL8)              :: beta_N
REAL(ZREAL8)              :: recip_alpha_f
REAL(ZREAL8)              :: recip_alpha_N
REAL(ZREAL8)              :: bdiva_f
REAL(ZREAL8)              :: bdiva_N
REAL(ZREAL8)              :: fourpi2
REAL(ZREAL8)              :: sr_nlongs
REAL(ZREAL8)              :: half_sr_nlongs




!**************************************************************************************************
! Declare Compound Types
!**************************************************************************************************

!-----------------------------------------------------------------------
! To store raw UM data slice
TYPE UM_type
  REAL(ZREAL8) :: longs_u(1:nlongs)
  REAL(ZREAL8) :: longs_v(1:nlongs)
  REAL(ZREAL8) :: half_levs(1:nlevs+1)
  REAL(ZREAL8) :: full_levs(0:nlevs)
  REAL(ZREAL8) :: u(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: v(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: w(1:nlongs,0:nlevs)
  REAL(ZREAL8) :: density(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: theta(1:nlongs,1:nlevs)
  REAL(ZREAL8) :: exner_pressure(1:nlongs,1:nlevs+1)
  REAL(ZREAL8) :: orog_height(1:nlongs)
END TYPE UM_type


!-----------------------------------------------------------------------
! To store information about the dimensions (axes)
TYPE dims_type
  REAL(ZREAL8) :: longs_u(0:nlongs+1)
  REAL(ZREAL8) :: longs_v(0:nlongs+1)
  REAL(ZREAL8) :: half_levs(0:nlevs+1)
  REAL(ZREAL8) :: full_levs(0:nlevs+1)
  ! Variables required for vertical interpolation
  REAL(ZREAL8) :: a1(1:nlevs+1)
  REAL(ZREAL8) :: b1(1:nlevs+1)
  REAL(ZREAL8) :: a2(0:nlevs)
  REAL(ZREAL8) :: b2(0:nlevs)
  REAL(ZREAL8) :: recip_half_kp1_k(0:nlevs)
  REAL(ZREAL8) :: recip_half_k_km1(1:nlevs+1)
  REAL(ZREAL8) :: recip_full_kp1_k(0:nlevs)
  REAL(ZREAL8) :: recip_full_k_km1(1:nlevs+1)
END TYPE dims_type


!-----------------------------------------------------------------------
! To store ABC model fields and some diagnostics (single time)
TYPE ABC_type
  ! Horizontal grid is Awakara C grid
  ! Vertical grid is Charney-Phillips
  REAL(ZREAL8) :: u(0:nlongs+1,0:nlevs+1)               ! Zonal wind perturbation
  REAL(ZREAL8) :: v(0:nlongs+1,0:nlevs+1)               ! Meridional wind perturbation
  REAL(ZREAL8) :: w(0:nlongs+1,0:nlevs+1)               ! Vertical wind perturbation
  REAL(ZREAL8) :: r(0:nlongs+1,0:nlevs+1)               ! rho density perturbation
  REAL(ZREAL8) :: b(0:nlongs+1,0:nlevs+1)               ! buoyancy perturbation
  REAL(ZREAL8) :: tracer(0:nlongs+1,0:nlevs+1)          ! Tracer
  REAL(ZREAL8) :: rho(0:nlongs+1,0:nlevs+1)             ! rho full field
  REAL(ZREAL8) :: b_ef(0:nlongs+1,0:nlevs+1)            ! Effective buoyancy
  REAL(ZREAL8) :: hydro_imbal(0:nlongs+1,0:nlevs+1)     ! Hydrostatic imbalance
  REAL(ZREAL8) :: geost_imbal(0:nlongs+1,0:nlevs+1)     ! Geostrophic imbalance
  REAL(ZREAL8) :: vert_mom_source(0:nlongs+1,0:nlevs+1) ! Vertical momentum source
  REAL(ZREAL8) :: horiz_div(0:nlongs+1,0:nlevs+1)       ! Horizontal divergence
  REAL(ZREAL8) :: horiz_vort(0:nlongs+1,0:nlevs+1)      ! Horizontal divergence
  REAL(ZREAL8) :: Kinetic_Energy
  REAL(ZREAL8) :: Buoyant_Energy
  REAL(ZREAL8) :: Elastic_Energy
  REAL(ZREAL8) :: Total_Energy
END TYPE ABC_type


!-----------------------------------------------------------------------
! Used with model integration scheme
TYPE Averages_type
  REAL(ZREAL8) :: u_1(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: u_2(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: u_m(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_1(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_2(0:nlongs+1, 0:nlevs+1)
  REAL(ZREAL8) :: w_m(0:nlongs+1, 0:nlevs+1)
END TYPE Averages_type


!-----------------------------------------------------------------------
! Scheme to store control variable fields (classic use commented)
TYPE CV_type
  REAL(ZREAL8) :: v1(0:nlongs+1,0:nlevs+1)  ! streamfunction
  REAL(ZREAL8) :: v2(0:nlongs+1,0:nlevs+1)  ! velocity potential
  REAL(ZREAL8) :: v3(0:nlongs+1,0:nlevs+1)  ! (unbalanced) r
  REAL(ZREAL8) :: v4(0:nlongs+1,0:nlevs+1)  ! (unbalanced) b
  REAL(ZREAL8) :: v5(0:nlongs+1,0:nlevs+1)  ! (unbalanced) w
  REAL(ZREAL8) :: v6(0:nlongs+1,0:nlevs+1)  ! tracer
END TYPE CV_type



!-----------------------------------------------------------------------
! Scheme to store the control variable transform data (CVT)
TYPE CVT_type

  ! Options for the transforms
  INTEGER      :: CVT_order          ! 1 = as original MetO
                                     ! 2 = reversed horiz/vert
                                     ! 3 = as REP's thesis
  INTEGER      :: CVT_param_opt_gb   ! 1 = analytical balance (geostrophic)
                                     ! 2 = statistical balance
                                     ! 3 = no geostrophic balance
  INTEGER      :: CVT_param_opt_hb   ! 1 = analytical balance (hydrostatic)
                                     ! 2 = statistical balance
                                     ! 3 = no hydrostatic balance
  INTEGER      :: CVT_param_opt_ab   ! 1 = analytical balance (anelastic for w)
                                     ! 2 = no anelastic balance
  INTEGER      :: CVT_param_opt_reg  ! 1 = use vertical regression of the gb r
                                     ! 2 = no vertical regression
  INTEGER      :: CVT_vert_opt_sym   ! 1 = non-symmetric transform
                                     ! 2 = symmetric transform

  INTEGER      :: CVT_stddev_opt     ! 1 = Stddev constant for each control variable
                                     ! 2 = Level dependent only
                                     ! 3 = Longitude and level dependent

  ! Data structures to hold the transforms
  ! Standard deviations of the 6 control parameters
  REAL(ZREAL8) :: sigma1(1:nlongs, 1:nlevs)
  REAL(ZREAL8) :: sigma2(1:nlongs, 1:nlevs)
  REAL(ZREAL8) :: sigma3(1:nlongs, 1:nlevs)
  REAL(ZREAL8) :: sigma4(1:nlongs, 1:nlevs)
  REAL(ZREAL8) :: sigma5(1:nlongs, 1:nlevs)
  REAL(ZREAL8) :: sigma6(1:nlongs, 1:nlevs)
  ! Vertical modes of the 6 control parameters
  REAL(ZREAL8) :: VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1)
  ! Vertical eigenvalues of the 6 control parameters (these are actually the square-roots)
  REAL(ZREAL8) :: VertEV1(1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertEV2(1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertEV3(1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertEV4(1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertEV5(1:nlevs, 1:nlongs/2+1)
  REAL(ZREAL8) :: VertEV6(1:nlevs, 1:nlongs/2+1)
  ! Horizontal eigenvalues of the 6 control parameters (these are actually the square-roots)
  REAL(ZREAL8) :: HorizEV1(1:nlongs/2+1, 1:nlevs)
  REAL(ZREAL8) :: HorizEV2(1:nlongs/2+1, 1:nlevs)
  REAL(ZREAL8) :: HorizEV3(1:nlongs/2+1, 1:nlevs)
  REAL(ZREAL8) :: HorizEV4(1:nlongs/2+1, 1:nlevs)
  REAL(ZREAL8) :: HorizEV5(1:nlongs/2+1, 1:nlevs)
  REAL(ZREAL8) :: HorizEV6(1:nlongs/2+1, 1:nlevs)
  ! Regression data for balanced density
  REAL(ZREAL8) :: Cov_rbalrbal(1:nlevs, 1:nlevs)
  REAL(ZREAL8) :: Cov_rtotrbal(1:nlevs, 1:nlevs)
  REAL(ZREAL8) :: Regression(1:nlevs, 1:nlevs)

END TYPE CVT_type


!-----------------------------------------------------------------------
! Storing the observations
TYPE Obs_type
  ! Batch index (for later development, e.g. correlated obs errs)
  INTEGER      :: batch
  INTEGER      :: obnumber_thisfile

  ! Observation time
  !INTEGER      :: year
  !INTEGER      :: month
  !INTEGER      :: day
  !INTEGER      :: hour
  !INTEGER      :: min
  !INTEGER      :: sec
  ! Absoute (seconds)
  INTEGER      :: t

  ! Location
  REAL(ZREAL8) :: longitude_deg
  REAL(ZREAL8) :: level_ht
  INTEGER      :: xbox_lower    ! Closest lower grid index (x)
  INTEGER      :: xbox_lower_ws ! Auxillary info needed for wind speed obs
  INTEGER      :: zbox_lower    ! Closest lower grid index (z)
  INTEGER      :: zbox_lower_ws ! Auxillary info needed for wind speed obs
  INTEGER      :: tstep_lower   ! Closest lower time index of corresponding model state sequence

  ! Observations
  LOGICAL      :: ob_ok
  INTEGER      :: ob_of_what    ! 1 (u), 2 (v), 3 (w), 4 (r), 5 (b), 6 (tracer)
                                ! 7 (horizontal wind speed)
                                ! 8 (total wind speed)
  LOGICAL      :: y_true_known
  REAL(ZREAL8) :: y_true        ! True value of ob
  REAL(ZREAL8) :: y             ! Observation value
  REAL(ZREAL8) :: stddev        ! Error standard deviation
  REAL(ZREAL8) :: variance      ! Error variance (stddev squared)
  REAL(ZREAL8) :: y_ref         ! Model observation (reference)
  REAL(ZREAL8) :: d             ! y - y_ref
  REAL(ZREAL8) :: deltay_m      ! Model observation increment
  REAL(ZREAL8) :: hxmy          ! deltay_m - d = y_ref + deltay_m - y
  REAL(ZREAL8) :: deltay_m_hat  ! d(JO)/d(deltay_m) = (R^-1) hxmy

  ! Pointer to the next observation record
  TYPE(Obs_type), POINTER :: next

END TYPE Obs_type


!-----------------------------------------------------------------------
! Specification of the observations for generation
TYPE ObsSpec_type
  ! Reference time for all observations
  INTEGER      :: year0
  INTEGER      :: month0
  INTEGER      :: day0
  INTEGER      :: hour0
  INTEGER      :: min0
  INTEGER      :: sec0

  INTEGER      :: NumBatches                  ! Number of obervation batches
  INTEGER      :: batch(1:maxbatches)         ! The batch numbers
  INTEGER      :: seconds(1:maxbatches)       ! The absolute time of this batch (seconds)
  INTEGER      :: ob_of_what(1:maxbatches)    ! What is observed (see same variable name in Obs_type)
  INTEGER      :: NumObs_long(1:maxbatches)   ! Number of observations in longitude direction
  INTEGER      :: NumObs_height(1:maxbatches) ! Number of observations in the height direction
  REAL(ZREAL8) :: long_min(1:maxbatches)      ! min longitude of obs patch
  REAL(ZREAL8) :: long_max(1:maxbatches)      ! max longitude of obs patch
  REAL(ZREAL8) :: height_min(1:maxbatches)    ! min height of obs patch
  REAL(ZREAL8) :: height_max(1:maxbatches)    ! max height of obs patch
  REAL(ZREAL8) :: stddev(1:maxbatches)        ! error standard dev
END TYPE ObsSpec_type





! Variables to do with preparing the ABC init state
! -------------------------------------------------
INTEGER                   :: Init_ABC_opt                 ! 1=Take UM data
                                                          ! 2=Pressure blob
CHARACTER(LEN=256)        :: datadirUM=''                 ! Directory to do with UM data
CHARACTER(LEN=256)        :: init_um_file=''              ! Input UM filename
INTEGER                   :: latitude = 144               ! The latitude to be extracted
LOGICAL                   :: Regular_vert_grid = .TRUE.   ! Set to use regular vertical levels
LOGICAL                   :: gravity_wave_switch = .FALSE.! To set u=0 (simulate gws)
CHARACTER(LEN=256)        :: init_ABC_file=''             ! Initial ABC filename
INTEGER                   :: press_source_x = 180         ! Long pos of centre of pressure blob
INTEGER                   :: press_source_z = 30          ! Level pos of centre of pressure blob
INTEGER                   :: x_scale = 80                 ! Long size of pressure blob
INTEGER                   :: z_scale = 3                  ! Level size of pressure blob

! Variables to do with running the forward model
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirABC_in=''             ! Directory to do with input of simplified model data
CHARACTER(LEN=256)        :: datadirABC_out=''            ! Directory to do with output of simplified model data
CHARACTER(LEN=256)        :: output_ABC_file=''           ! Dump file
CHARACTER(LEN=256)        :: diagnostics_file=''          ! For diagnostics
REAL(ZREAL8)              :: runlength = 60.0             ! Runlength in seconds
INTEGER                   :: ndumps = 10                  ! number of dump times
LOGICAL                   :: convection_switch = .FALSE.  ! Set to ?
REAL(ZREAL8)              :: press_amp = 0.01             ! Amplitude of pressure blob
LOGICAL                   :: Adv_tracer = .FALSE.         ! Set to advect tracer in calculations
LOGICAL                   :: Lengthscale_diagnostics = .FALSE.
                                                          ! Set to do lengthscale diagnostics at final time

! Variables to do with the calibration of the CVT, etc.
! ----------------------------------------
INTEGER                   :: CalibRunStage = 1            ! Stage 1 is to convert UM to ABC forecasts
                                                          ! Stage 2 is to compute perturbations
                                                          ! Stage 3 is to determine the regression
                                                          ! Stage 4 is parameter transform
                                                          ! Stage 5 is calibration of spatial stats
INTEGER                   :: NEns = 50                    ! Number of ensembles (0=do not use ensembles)
CHARACTER(LEN=256)        :: EnsDirs(1:Nensmax)=''        ! Directories containing the ensembles
INTEGER                   :: NEnsMems = 24                ! Number of ensemble members

INTEGER                   :: NNMC = 0                     ! Number of NMC pairs (0=do not use NMC)
CHARACTER(LEN=256)        :: NMCDirs(1:NNMCmax)=''        ! Directories containing the NMC pairs

INTEGER                   :: Nlats = 1                    ! Number of different lats to extract from
                                                          ! UM files
INTEGER                   :: latindex(1:Nlatsmax)         ! The latitude indices

CHARACTER(LEN=256)        :: datadirABCfcs=''             ! Directory of ABC model forecasts (stage 1 calibration)
CHARACTER(LEN=256)        :: datadirABCperts=''           ! Directory of ABC model forecasts (stage 2 calibration)
CHARACTER(LEN=256)        :: datadirRegression=''         ! Directory of regression diagnostics (stage 3 calibration)
CHARACTER(LEN=256)        :: datadirConParams=''          ! Directory of control parameters (stage 4 calibration)
CHARACTER(LEN=256)        :: datadirCVT=''                ! Directory of the CVT data (stage 5 calibration)
CHARACTER(LEN=256)        :: CVT_file=''                  ! Name of file containing CVT data

TYPE(CVT_type)            :: CVT                          ! Default values are set inside SetOptions

INTEGER                   :: VertSmoothPoints = 0         ! Number of points in vertical to average for standard dev.
INTEGER                   :: HorizSmoothPoints = 0        ! Number of points in horizontal to average for standard dev.


! Variables to do with testing the DA
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirTestDA=''             ! Directory to do with testing the DA components
CHARACTER(LEN=256)        :: LS_file=''                   ! Name of file containing LS
CHARACTER(LEN=256)        :: Pert_file=''                 ! Name of file containing pert data


LOGICAL                   :: RunAdjTests_CVT = .FALSE.
LOGICAL                   :: RunAdjTests_obs = .FALSE.
LOGICAL                   :: RunInvTests     = .FALSE.


! Variables to do with generating a background state, and the observations
! ----------------------------------------------
INTEGER                   :: Generate_mode = 1            ! 1 = Generate file that specifies obs positions/times/etc
                                                          ! 2 = Generate synthetic observations
                                                          ! 3 = Generate synthetic background state
TYPE(ObsSpec_type)        :: ObsSpec                      ! Initialised in SetOptions
CHARACTER(LEN=256)        :: datadir_ObsSpec=''           ! Directory to place obs spec file
CHARACTER(LEN=256)        :: ObsSpec_file=''              ! Filename of observation specifications
CHARACTER(LEN=256)        :: datadir_Bg=''                ! Directory containing background
CHARACTER(LEN=256)        :: Bg_file=''                   ! Background state filename
CHARACTER(LEN=256)        :: datadir_Obs=''               ! Directory containing observations
CHARACTER(LEN=256)        :: Obs_file=''                  ! Observations filename
REAL(ZREAL8)              :: dt_da = 60.0                 ! Time-step of the DA (seconds)


! Variables to do with implied/raw covariances
! ----------------------------------------------
INTEGER                   :: ImplCov_npoints = 0          ! Number of different source points to consider
INTEGER                   :: longindex(1:Npointsmax)      ! The longitude indices
INTEGER                   :: levindex(1:Npointsmax)       ! The level indices
CHARACTER(LEN=256)        :: datadirImpliedCov=''         ! Directory to output implied covariances
CHARACTER(LEN=256)        :: datadirRawCov=''             ! Directory to output raw covariances


! Variables to do with the DA
! ----------------------------------------------
CHARACTER(LEN=256)        :: datadirAnal=''               ! Directory for the analysis files
CHARACTER(LEN=256)        :: anal_file=''                 ! Analysis file
CHARACTER(LEN=256)        :: analinc_file=''              ! Analysis increment file
INTEGER                   :: t0 = 0                       ! Time of start of this DA cycle (seconds)
INTEGER                   :: Hybrid_opt = 1               ! 1 = standard B
                                                          ! 2 = pure EnVar
                                                          ! 3 = hybrid EnVar
                                                          ! 4 = reduced rank KF-type hybrid
INTEGER                   :: Vartype = 3                  ! 3 = 3DVar
                                                          ! 35 = 3DFGAT
                                                          ! 4 = 4DVar
INTEGER                   :: N_outerloops = 1             ! Number of outer loops
INTEGER                   :: N_innerloops_max = 10        ! Maximum number of inner loops
REAL(ZREAL8)              :: mu = 0.001                   ! Small number for perturbing search direction
REAL(ZREAL8)              :: minus_mu                     ! Negative of above
REAL(ZREAL8)              :: crit_inner = 0.01            ! Stopping criterion for inner loop

! Variables to do with FFTs
! -------------------------
INTEGER, PARAMETER        :: fft_worklen_x = 2*nlongs     ! Much larger than needbe
                                                          ! At least nlongs + INT(LOG(REAL(nlongs,kind=8)) / LOG(2.0)) + 5
LOGICAL                   :: fft_init_x = .FALSE.
REAL(ZREAL8)              :: fft_wsave_x(1:fft_worklen_x)
REAL(ZREAL8)              :: fft_work_x(1:nlongs)
INTEGER, PARAMETER        :: fft_worklen_z = 2*nlevs      ! Much longer than needbe
                                                          ! At least INT(LOG(REAL(nlevs,kind=8)) / LOG(2.0)) + 5
LOGICAL                   :: fft_init_z = .FALSE.
REAL(ZREAL8)              :: fft_wsave_z(1:fft_worklen_z)
REAL(ZREAL8)              :: fft_work_z(1:nlevs)


! Variables to do with linear analysis of the model
! -------------------------------------------------
CHARACTER(LEN=256)        :: datadirLinearAnal=''         ! Separate directory for the linear analysis data





! Define namelist
! ---------------
NAMELIST / UserOptions /                                                           &
! --- Set initial state for ABC model ---
  Init_ABC_opt, datadirUM, init_um_file, latitude, Regular_vert_grid,              &
  gravity_wave_switch, f, A, B, C, BoundSpread, init_ABC_file,                     &
! --- Run the forward model ---
  datadirABC_in, datadirABC_out, output_ABC_file, diagnostics_file, dt, dx, H,     &
  runlength, ndumps,                                                               &
  convection_switch, press_source_x, press_source_z, x_scale, z_scale, press_amp,  &
  Adv_tracer, Lengthscale_diagnostics,                                             &
! --- Testing the DA ---
  datadirTestDA, RunAdjTests_CVT, RunAdjTests_obs, RunInvTests, LS_file, Pert_file,&
! --- DA ---
  datadirCVT, CVT_file, Hybrid_opt, Vartype, datadirAnal, anal_file, analinc_file, &
  N_outerloops, N_innerloops_max, crit_inner,                                      &
! --- Linear analysis ---
  datadirLinearAnal,                                                               &
! --- Calibration, CVT, etc.
  CalibRunStage, NEns, EnsDirs, NEnsMems, NNMC, NMCDirs, datadirConParams,         &
  datadirABCfcs, datadirABCperts, datadirRegression,                               &
  Nlats, latindex, CVT, VertSmoothPoints, HorizSmoothPoints,                       &
! --- Generate background, obs, etc.
  Generate_mode, ObsSpec, datadir_ObsSpec, ObsSpec_file, datadir_Bg,               &
  datadir_Obs, Bg_file, Obs_file, dt_da, t0, random_seed,                          &
! Implied covs
  datadirImpliedCov, datadirRawCov, ImplCov_npoints, longindex, levindex



END MODULE DefConsTypes
