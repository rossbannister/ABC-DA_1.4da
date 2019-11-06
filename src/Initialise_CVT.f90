SUBROUTINE Initialise_CVT (CVT)

! Initialise control variable transform structure
! Does not alter the ABC parameters and options though

USE DefConsTypes, ONLY :     &
    CVT_type,                &
    nlongs, nlevs

IMPLICIT NONE

TYPE(CVT_type), INTENT(INOUT)   :: CVT

  ! Standard deviations of the 6 control parameters
  CVT % sigma1(1:nlongs, 1:nlevs)                 = 0.0
  CVT % sigma2(1:nlongs, 1:nlevs)                 = 0.0
  CVT % sigma3(1:nlongs, 1:nlevs)                 = 0.0
  CVT % sigma4(1:nlongs, 1:nlevs)                 = 0.0
  CVT % sigma5(1:nlongs, 1:nlevs)                 = 0.0
  CVT % sigma6(1:nlongs, 1:nlevs)                 = 0.0

  ! Vertical modes of the 6 control parameters
  CVT % VertMode1(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode2(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode3(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode4(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode5(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0
  CVT % VertMode6(1:nlevs, 1:nlevs, 1:nlongs/2+1) = 0.0

  ! Vertical eigenvalues of the 6 control parameters (these are actually the square-roots)
  CVT % VertEV1(1:nlevs, 1:nlongs/2+1)            = 0.0
  CVT % VertEV2(1:nlevs, 1:nlongs/2+1)            = 0.0
  CVT % VertEV3(1:nlevs, 1:nlongs/2+1)            = 0.0
  CVT % VertEV4(1:nlevs, 1:nlongs/2+1)            = 0.0
  CVT % VertEV5(1:nlevs, 1:nlongs/2+1)            = 0.0
  CVT % VertEV6(1:nlevs, 1:nlongs/2+1)            = 0.0

  ! Horizontal eigenvalues of the 6 control parameters (these are actually the square-roots)
  CVT % HorizEV1(1:nlongs/2+1, 1:nlevs)           = 0.0
  CVT % HorizEV2(1:nlongs/2+1, 1:nlevs)           = 0.0
  CVT % HorizEV3(1:nlongs/2+1, 1:nlevs)           = 0.0
  CVT % HorizEV4(1:nlongs/2+1, 1:nlevs)           = 0.0
  CVT % HorizEV5(1:nlongs/2+1, 1:nlevs)           = 0.0
  CVT % HorizEV6(1:nlongs/2+1, 1:nlevs)           = 0.0

  ! Regression data for balanced density
  CVT % Cov_rbalrbal(1:nlevs, 1:nlevs)            = 0.0
  CVT % Cov_rtotrbal(1:nlevs, 1:nlevs)            = 0.0
  CVT % Regression(1:nlevs, 1:nlevs)              = 0.0

END SUBROUTINE Initialise_CVT
