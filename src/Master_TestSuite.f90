PROGRAM Master_TestSuite

!*****************************************************
!*   Code to do the following tests for the DA       *
!*                                                   *
!*   Adjoint tests                                   *
!*   Inverse tests                                   *
!*   Linearisation tests (not yet implemented)       *
!*                                                   *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    ZREAL8,                      &
    nlongs, nlevs,               &
    dims_type,                   &
    ABC_type,                    &
    CV_type,                     &
    CVT_type,                    &
    datadirTestDA,               &
    datadirABCfcs,               &
    datadirABCperts,             &
    datadirCVT,                  &
    LS_file,                     &
    Pert_file,                   &
    CVT_file,                    &
    RunAdjTests_CVT,             &
    RunAdjTests_obs,             &
    RunInvTests,                 &
    diagnostics_file,            &
    datadir_Obs,                 &
    Obs_file,                    &
    Obs_type



IMPLICIT NONE

INCLUDE "Boundaries.interface"
INCLUDE "Boundaries_adj.interface"
INCLUDE "InnerProdModelSpace.interface"
INCLUDE "InnerProdControlSpace.interface"
INCLUDE "Boundaries_CV.interface"
INCLUDE "Boundaries_CV_adj.interface"
INCLUDE "Read_Obs.interface"
INCLUDE "ModelObservations.interface"
INCLUDE "ModelObservations_linear.interface"
INCLUDE "InnerProduct_obsself.interface"
INCLUDE "ModelObservations_adj.interface"
INCLUDE "Write_Obs.interface"
INCLUDE "DeAllocate_Obs.interface"


! Declare variables
!==========================
TYPE(dims_type)             :: dims
TYPE(ABC_type)              :: Model_data1, Model_data2, LS, Pert1, Pert2, Pert3
TYPE(CV_type)               :: CV_data1, CV_data2, CV_data3
TYPE(ABC_type)              :: Model_data_array1(1:2), Model_data_array2(1:2)
TYPE(ABC_type)              :: LSarray(1:2)
TYPE(CVT_type)              :: CVT
REAL(ZREAL8)                :: Field1(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: Field1a(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: Field2(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: Field2a(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: Field3(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: Field3a(0:nlongs+1,0:nlevs+1)
REAL(ZREAL8)                :: TwoPoints1(1:2), TwoPoints2(1:2)
REAL(ZREAL8)                :: AxisValuesx(1:2), AxisValuest(1:2), ObValue
REAL(ZREAL8)                :: Value
REAL(ZREAL8)                :: LHS, RHS
REAL(ZREAL8)                :: testmatrix(1:nlevs,1:nlevs), check(1:nlevs,1:nlevs)
COMPLEX(ZREAL8)             :: LHScomplex, RHScomplex
LOGICAL                     :: switches(1:6)
INTEGER                     :: loop, z

CHARACTER(LEN=320)          :: LS_filename, Pert_filename, Test_diags_filename, CVT_filename
CHARACTER(LEN=320)          :: InvTestOutput_filename, Obs_filename, ABC_init_filename
TYPE(Obs_type), POINTER     :: Observations
REAL(ZREAL8)                :: dt, dt_da
INTEGER                     :: timestepsreq, DAtimestepsreq, maxtime
TYPE(ABC_type), ALLOCATABLE :: LSfc(:), ABC1(:), ABC2(:)

REAL(ZREAL8)                :: INT_HF, INT_FH
REAL(ZREAL8)                :: Interpolate1D, Interpolate3D


PRINT*, '*************************************************************************'
PRINT*, 'Running Master_TestSuite'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions


  LS_filename          = TRIM(datadirABCfcs)   // '/' // TRIM(LS_file)
  Pert_filename        = TRIM(datadirABCperts) // '/' // TRIM(Pert_file)
  CVT_filename         = TRIM(datadirCVT)      // '/' // TRIM(CVT_file)
  Test_diags_filename  = TRIM(datadirTestDA)   // '/' // TRIM(diagnostics_file)


  OPEN (12, file=Test_diags_filename)

  ! Read in LS
  PRINT*, 'Reading in linearisation state ...'
  CALL Read_state_2d (LS_filename, LS, dims, 1)
  PRINT*, '-- done'

  ! Set some commonly-used constants
  CALL Set_ht_dep_cons (dims)


  IF (RunAdjTests_CVT .OR. RunInvTests) THEN
    ! Read in CVT data prepared beforehand
    PRINT*, 'Reading in CVT data ...'
    CVT_filename = TRIM(datadirCVT) // '/' // TRIM(CVT_file)
    CALL Read_Covs (CVT_filename, CVT,              &
                    .TRUE., .TRUE., .TRUE., .TRUE., .TRUE.)
    PRINT *, '-- done'
  END IF



  IF (RunAdjTests_CVT) THEN
    PRINT*, '===== Running adjoint tests for CVT ====='
    WRITE (12,*) '===== Running adjoint tests for CVT ====='
    WRITE (12,*) '(A v)^T A v  =?  (A^T A v)^T v'
    ! --------------------------------------------------------------------------


    PRINT *, 'Boundaries'
    WRITE (12,*) '--- Boundaries ---'
    DO loop = 1, 7
      IF (loop < 7) THEN
        switches(1:6) = .FALSE.
        switches(loop) = .TRUE.
      ELSE
        switches(1:6) = .TRUE.
      END IF

      WRITE (12,*) 'Variables considered : ', switches(1:6)
      ! Put random numbers in a model data structure
      CALL Initialise_model_vars (Model_data1, .TRUE.)                                                ! v

      ! Act with forward operator
      Model_data2 = Model_data1
      CALL Boundaries (Model_data2, set_u=switches(1), set_v=switches(2), set_w=switches(3),      &
                                    set_r=switches(4), set_b=switches(5), set_tracer=switches(6))     ! A v
      ! Compute the LHS
      LHS = InnerProdModelSpace (Model_data2, Model_data2,                                        &
                                 do_u=switches(1), do_v=switches(2), do_w=switches(3),            &
                                 do_r=switches(4), do_b=switches(5), do_tracer=switches(6))           ! (A v)^T A v

      ! Act with adjoint operator
      CALL Boundaries_adj (Model_data2, set_u=switches(1), set_v=switches(2), set_w=switches(3),  &
                                        set_r=switches(4), set_b=switches(5), set_tracer=switches(6)) ! A^T A v
      ! Compute the RHS
      RHS = InnerProdModelSpace (Model_data1, Model_data2,                                        &
                                 do_u=switches(1), do_v=switches(2), do_w=switches(3),            &
                                 do_r=switches(4), do_b=switches(5), do_tracer=switches(6))           ! (A^T A v)^T v

      WRITE (12,*) 'LHS = ', LHS
      WRITE (12,*) 'RHS = ', RHS
    END DO


    ! --------------------------------------------------------------------------
    PRINT *, 'Boundaries_CV'
    WRITE (12,*) '--- Boundaries_CV ---'
    DO loop = 1, 7
      IF (loop < 7) THEN
        switches(1:6) = .FALSE.
        switches(loop) = .TRUE.
      ELSE
        switches(1:6) = .TRUE.
      END IF

      WRITE (12,*) 'Variables considered : ', switches(1:6)
      ! Put random numbers in a model data structure
      CALL Initialise_CVs (CV_data1, .TRUE.)                                                     ! v

      ! Act with forward operator
      CV_data2 = CV_data1
      CALL Boundaries_CV (CV_data2, set_1=switches(1), set_2=switches(2), set_3=switches(3),      &
                                    set_4=switches(4), set_5=switches(5), set_6=switches(6))     ! A v
      ! Compute the LHS
      LHScomplex = InnerProdControlSpace (CV_data2, CV_data2,                                     &
                                   do_1=switches(1), do_2=switches(2), do_3=switches(3),          &
                                   do_4=switches(4), do_5=switches(5), do_6=switches(6))         ! (A v)^T A v
      ! Act with adjoint operator
      CALL Boundaries_CV_adj (CV_data2, set_1=switches(1), set_2=switches(2), set_3=switches(3),  &
                                        set_4=switches(4), set_5=switches(5), set_6=switches(6)) ! A^T A v
      ! Compute the RHS
      RHScomplex = InnerProdControlSpace (CV_data1, CV_data2,                                     &
                                   do_1=switches(1), do_2=switches(2), do_3=switches(3),          &
                                   do_4=switches(4), do_5=switches(5), do_6=switches(6))         ! (A^T A v)^T v

      WRITE (12,*) 'LHS = ', LHScomplex
      WRITE (12,*) 'RHS = ', RHScomplex
    END DO


    ! --------------------------------------------------------------------------
    PRINT *, 'LinearBal_r'
    WRITE (12,*) '--- LinearBal_r ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(0:nlongs+1,1:nlevs))                             ! v

    ! Act with forward operator
    CALL LinearBal_r (Field1(0:nlongs+1,1:nlevs), Field2(1:nlongs,1:nlevs))     ! A v
    ! Compute the LHS
    LHS = SUM(Field2(1:nlongs,1:nlevs) * Field2(1:nlongs,1:nlevs))              !(A v)^T A v

    ! Act with adjoint operator
    Field3(0:nlongs+1,1:nlevs) = 0.0
    CALL LinearBal_r_adj (Field3(0:nlongs+1,1:nlevs), Field2(1:nlongs,1:nlevs)) ! A^T A v
    ! Compute the RHS
    RHS = SUM(Field3(0:nlongs+1,1:nlevs) * Field1(0:nlongs+1,1:nlevs))          ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'Anbalw'
    WRITE (12,*) '--- Anbalw ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(0:nlongs+1,0:nlevs+1))                            ! v

    ! Act with forward operator
    CALL Anbalw (LS % rho(0:nlongs+1,0:nlevs+1), Field1(0:nlongs+1,0:nlevs+1),      &
                 Field2(1:nlongs,1:nlevs), dims)                                 ! A v
    ! Compute the LHS
    LHS = SUM(Field2(1:nlongs,1:nlevs) * Field2(1:nlongs,1:nlevs))               ! (A v)^T (A v)

    ! Act with adjoint operator
    Field3(0:nlongs+1,0:nlevs+1) = 0.0
    CALL Anbalw_adj (LS % rho(0:nlongs+1,0:nlevs+1), Field3(0:nlongs+1,0:nlevs+1),  &
                     Field2(1:nlongs,1:nlevs), dims)                             ! A^T A v
    ! Compute the RHS
    RHS = SUM(Field3(0:nlongs+1,0:nlevs+1) * Field1(0:nlongs+1,0:nlevs+1))       ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'Helmholtz'
    WRITE (12,*) '--- Helmholtz ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(0:nlongs+1,1:nlevs))                            ! v
    CALL RANDOM_NUMBER (Field1a(0:nlongs+1,1:nlevs))                           ! v

    ! Act with forward operator
    CALL Helmholtz (Field1(0:nlongs+1,1:nlevs), Field1a(0:nlongs+1,1:nlevs),     &
                    Field2(1:nlongs,1:nlevs), Field2a(1:nlongs,1:nlevs))       ! A v

    ! Compute the LHS
    LHS = SUM(Field2(1:nlongs,1:nlevs) * Field2(1:nlongs,1:nlevs)) +             &
          SUM(Field2a(1:nlongs,1:nlevs) * Field2a(1:nlongs,1:nlevs))           ! (A v)^T (A v)

    ! Act with adjoint operator
    Field3(0:nlongs+1,1:nlevs)  = 0.0
    Field3a(0:nlongs+1,1:nlevs) = 0.0
    CALL Helmholtz_adj (Field3(0:nlongs+1,1:nlevs), Field3a(0:nlongs+1,1:nlevs), &
                        Field2(1:nlongs,1:nlevs), Field2a(1:nlongs,1:nlevs))   ! A^T A v
    ! Compute the RHS
    RHS = SUM(Field3(0:nlongs+1,1:nlevs) * Field1(0:nlongs+1,1:nlevs)) +         &
          SUM(Field3a(0:nlongs+1,1:nlevs) * Field1a(0:nlongs+1,1:nlevs))       ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'INT_HF'
    WRITE (12,*) '--- INT_HF ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(0:1,0))                                 ! v

    ! Act with forward operator
    Field2(0,0) = INT_HF (Field1(0,0), Field1(1,0), 6, dims)           ! A v

    ! Compute the LHS
    LHS = Field2(0,0) * Field2(0,0)                                    ! (A v)^T (A v)

    ! Act with adjoint operator
    Field3(0:1,0)  = 0.0
    CALL INT_HF_adj (Field3(0,0), Field3(1,0), Field2(0,0), 6, dims)   ! A^T A v

    ! Compute the RHS
    RHS = Field3(0,0) * Field1(0,0) + Field3(1,0) * Field1(1,0)        ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'INT_FH'
    WRITE (12,*) '--- INT_FH ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(0:1,0))                                 ! v

    ! Act with forward operator
    Field2(0,0) = INT_FH (Field1(0,0), Field1(1,0), 6, dims)           ! A v

    ! Compute the LHS
    LHS = Field2(0,0) * Field2(0,0)                                    ! (A v)^T (A v)

    ! Act with adjoint operator
    Field3(0:1,0)  = 0.0
    CALL INT_FH_adj (Field3(0,0), Field3(1,0), Field2(0,0), 6, dims)   ! A^T A v

    ! Compute the RHS
    RHS = Field3(0,0) * Field1(0,0) + Field3(1,0) * Field1(1,0)        ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'HydroBal_b'
    WRITE (12,*) '--- HydroBal_b ---'
    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (Field1(1:nlongs,0:nlevs+1))                                   ! v

    ! Act with forward operator
    CALL HydroBal_b (Field1(1:nlongs,0:nlevs+1), Field2(1:nlongs,1:nlevs), dims)      ! A v

    ! Compute the LHS
    LHS = SUM(Field2(1:nlongs,1:nlevs) * Field2(1:nlongs,1:nlevs))                    ! (A v)^T (A v)

    ! Act with adjoint operator
    Field3(1:nlongs,0:nlevs+1)  = 0.0
    CALL HydroBal_b_adj (Field3(1:nlongs,0:nlevs+1), Field2(1:nlongs,1:nlevs), dims)  ! A^T A v

    ! Compute the RHS
    RHS = SUM(Field3(1:nlongs,0:nlevs+1) * Field1(1:nlongs,0:nlevs+1))                ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'U_p'
    WRITE (12,*) '--- U_p ---'
    ! Put random numbers in a model data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                                 ! v
    CALL Boundaries_CV (CV_data1)

    ! Act with forward operator
    CALL U_p (LS, CV_data1, Model_Data1, CVT % CVT_order,                                 &
              CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,     &
              CVT % CVT_param_opt_reg,                                                    &
              CVT % Regression(1:nlevs,1:nlevs), dims)                     ! A v
    ! Compute the LHS
    LHS = InnerProdModelSpace (Model_data1, Model_data1,                                  &
                               do_u=.TRUE., do_v=.TRUE., do_w=.TRUE.,                     &
                               do_r=.TRUE., do_b=.TRUE., do_tracer=.TRUE.) ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_CVs (CV_data2, .FALSE.)
    CALL U_p_adj (LS, CV_data2, Model_Data1, CVT % CVT_order,                             &
                  CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab, &
                  CVT % CVT_param_opt_reg,                                                &
                  CVT % Regression(1:nlevs,1:nlevs), dims)                 ! A^T A v
    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data2, CV_data1,                               &
                                 do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                   &
                                 do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.)    ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHScomplex



    ! --------------------------------------------------------------------------
    PRINT *, 'U_v'
    WRITE (12,*) '--- U_v ---'
    ! Put random numbers in a model data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                         ! v

    ! Act with forward operator
    CALL U_v (CV_data1, CV_data2, CVT)                             ! A v
    ! Compute the LHS
    LHScomplex = InnerProdControlSpace (CV_data2, CV_data2,                        &
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                 &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.,                 &
                            ignore_halos=.TRUE.)                   ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_CVs (CV_data3, .FALSE.)
    CALL U_v_adj (CV_data3, CV_data2, CVT)                         ! A^T A v
    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data3, CV_data1,                        &
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                 &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.,                 &
                            ignore_halos=.TRUE.)                   !(A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHScomplex
    WRITE (12,*) 'RHS = ', RHScomplex



    ! --------------------------------------------------------------------------
    PRINT *, 'U_stddev'
    WRITE (12,*) '--- U_stddev ---'
    ! Put random numbers in a model data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                         ! v

    ! Act with forward operator
    CV_data2 = CV_data1
    CALL U_stddev (CV_data2, CVT)                                  ! A v
    ! Compute the LHS
    LHScomplex = InnerProdControlSpace (CV_data2, CV_data2,                        &
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                 &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.) ! (A v)^T A v

    ! Act with adjoint operator
    CV_data3 = CV_data2
    CALL U_stddev (CV_data3, CVT)                                  ! A^T A v
    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data3, CV_data1,                        &
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                 &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.) !(A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHScomplex
    WRITE (12,*) 'RHS = ', RHScomplex



    ! --------------------------------------------------------------------------
    PRINT *, 'FFT (real 2 spectral)'
    WRITE (12,*) '--- FFT (real to spectral) ---'
    CALL Initialise_CVs (CV_data1, .FALSE.)
    ! Put random numbers in a model data structure
    CALL Initialise_model_vars (Pert2, .TRUE.)                      ! v

    ! Act with forward operator
    CALL fft_real2spec (Pert2 % u(1:nlongs,1:nlevs), CV_data1 % v1(1:nlongs,1:nlevs))

    ! Compute the LHS
    LHScomplex = InnerProdControlSpace (CV_data1, CV_data1,                        &
                                        ComplexSpace=.TRUE.,                       &
                                        do_1=.TRUE.,                               &
                                        ignore_halos=.TRUE.)        ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_model_vars (Pert3, .FALSE.)
    CALL fft_spec2real (CV_data1 % v1(1:nlongs,1:nlevs), Pert3 % u(1:nlongs,1:nlevs))

    ! Compute the RHS
    RHS        = InnerProdModelSpace (Pert3, Pert2,                                &
                                      do_u=.TRUE.,                                 &
                                      ignore_halos=.TRUE.)          ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHScomplex
    WRITE (12,*) 'RHS = ', RHS


    ! --------------------------------------------------------------------------
    PRINT *, 'FFT (spectral to real)'
    WRITE (12,*) '--- FFT (spectral to real) ---'
    CALL Initialise_CVs (CV_data2, .FALSE.)
    ! Put random numbers in a CV data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                          ! v

    ! Act with forward operator
    CALL fft_spec2real (CV_data1 % v1(1:nlongs,1:nlevs), CV_data2 % v1(1:nlongs,1:nlevs))

    ! Compute the LHS
    LHScomplex = InnerProdControlSpace (CV_data2, CV_data2,                        &
                                        ComplexSpace=.FALSE.,                      &
                                        do_1=.TRUE.,                               &
                                        ignore_halos=.TRUE.)        ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_CVs (CV_data3, .FALSE.)
    CALL fft_real2spec (CV_data2 % v1(1:nlongs,1:nlevs), CV_data3 % v1(1:nlongs,1:nlevs))

    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data3, CV_data1,                        &
                                        ComplexSpace=.TRUE.,                       &
                                        do_1=.TRUE.,                               &
                                        ignore_halos=.TRUE.)        ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHScomplex
    WRITE (12,*) 'RHS = ', RHScomplex



    ! --------------------------------------------------------------------------
    PRINT *, 'U_h'
    WRITE (12,*) '--- U_h ---'
    CALL Initialise_CVs (CV_data2, .FALSE.)
    ! Put random numbers in a CV data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                          ! v

    ! Act with forward operator
    CALL U_h (CV_data1, CV_data2, CVT)                              ! A v
    ! Compute the LHS
    LHScomplex = InnerProdControlSpace (CV_data2, CV_data2,                        &
                        do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                     &
                        do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.,                     &
                        ignore_halos=.TRUE.)                        ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_CVs (CV_data3, .FALSE.)
    CALL U_h_adj (CV_data3, CV_data2, CVT)                          ! A^T A v
    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data3, CV_data1, ComplexSpace=.TRUE.,   &
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,                 &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.,                 &
                            ignore_halos=.TRUE.)                    ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHScomplex
    WRITE (12,*) 'RHS = ', RHScomplex


    ! --------------------------------------------------------------------------
    PRINT *, 'U_trans'
    WRITE (12,*) '--- U_trans ---'
      CALL Initialise_model_vars (Model_data1, .FALSE.)
    ! Put random numbers in a CV data structure
    CALL Initialise_CVs (CV_data1, .TRUE.)                                 ! v

    ! Act with forward operator
    CALL U_trans (LS, CV_data1, Model_data1, CVT, dims)                    ! A v
    ! Compute the LHS
    LHS = InnerProdModelSpace (Model_data1, Model_data1,                     &
                               do_u=.TRUE., do_v=.TRUE., do_w=.TRUE.,        &
                               do_r=.TRUE., do_b=.TRUE., do_tracer=.TRUE.) ! (A v)^T A v

    ! Act with adjoint operator
    CALL Initialise_CVs (CV_data2, .FALSE.)
    CALL U_trans_adj (LS, CV_data2, Model_data1, CVT, dims)                ! A^T A v
    ! Compute the RHS
    RHScomplex = InnerProdControlSpace (CV_data2, CV_data1, ComplexSpace=.TRUE.,&
                            do_1=.TRUE., do_2=.TRUE., do_3=.TRUE.,           &
                            do_4=.TRUE., do_5=.TRUE., do_6=.TRUE.)         ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHScomplex

  END IF



! =============================================================================

  IF (RunAdjTests_obs) THEN
    PRINT*, '===== Running adjoint tests for obs operators ====='
    WRITE (12,*) '===== Running adjoint tests for obs operators ====='
    WRITE (12,*) '(A v)^T A v  =?  (A^T A v)^T v'
    ! --------------------------------------------------------------------------


    ! --------------------------------------------------------------------------
    PRINT *, 'Interpolation in 1d'
    WRITE (12,*) '--- Interpolation in 1d ---'
    ! Set-up the axis
    AxisValuesx(1) = 6.5
    AxisValuesx(2) = 8.7
    Value          = 8.1

    ! Put random numbers in a model data structure
    CALL RANDOM_NUMBER (TwoPoints1(1:2))                                     ! v

    ! Act with forward operator
    ObValue  = Interpolate1D (TwoPoints1(1:2), AxisValuesx(1:2), Value)      ! Av
    ! Compute the LHS
    LHS      = ObValue * ObValue                                             ! (A v)^T A v
    ! Act with adjoint operator
    TwoPoints2(1:2) = 0.0
    CALL Interpolate1D_adj (ObValue, TwoPoints2(1:2), AxisValuesx(1:2), Value) ! A^T A v
    ! Compute the RHS
    RHS      = TwoPoints2(1) * TwoPoints1(1) + TwoPoints2(2) * TwoPoints1(2) ! (A^T A v)^T v

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS



    ! --------------------------------------------------------------------------
    ! Interpolation in 3d (1 time step)

    ! Set-up the axis for time
    AxisValuest(1) = 18.0

    ! Put random numbers in a model data structure
    CALL Initialise_model_vars (Model_data_array1(1), .TRUE.)                                                ! v

    DO loop = 1, 6
      PRINT *, 'Interpolation in 3d (1 time step, quantity', loop, ')'
      WRITE (12,*) '--- Interpolation in 3d (1 time step, quantity', loop, ') ---'

      ! Act with forward operator
      ObValue  = Interpolate3D (1, Model_data_array1(1:1),          &
                                dims % longs_u(4:5),                &
                                dims % half_levs(8:9),              &
                                AxisValuest(1:1),                   &
                                4, 8,                               &
                                dims % longs_u(4) + 0.004,          &
                                dims % half_levs(8) + 45.0,         &
                                AxisValuest(1),                     &
                                loop )                                     ! Av
      ! Compute the LHS
      LHS      = ObValue * ObValue                                         ! (A v)^T A v

      ! Act with adjoint operator
      ! Initialize model data structure
      CALL Initialise_model_vars (Model_data_array2(1), .FALSE.)
      CALL Interpolate3D_adj (ObValue, 1, Model_data_array2(1:1),   &
                              dims % longs_u(4:5),                  &
                              dims % half_levs(8:9),                &
                              AxisValuest(1:1),                     &
                              4, 8,                                 &
                              dims % longs_u(4) + 0.004,            &
                              dims % half_levs(8) + 45.0,           &
                              AxisValuest(1),                       &
                              loop )                                       ! A^T A v
      ! Compute the RHS
      RHS      = InnerProdModelSpace (Model_data_array1(1),         &
                                      Model_data_array2(1),         &
                                      ignore_halos=.TRUE.)                 ! (A^T A v)^T v
      WRITE (12,*) 'LHS = ', LHS
      WRITE (12,*) 'RHS = ', RHS
    END DO



    ! --------------------------------------------------------------------------
    ! Interpolation in 3d (2 time steps)

    ! Set-up the axis for time
    AxisValuest(1) = 18.0
    AxisValuest(2) = 20.5

    ! Put random numbers in a model data structure
    CALL Initialise_model_vars (Model_data_array1(1), .TRUE.)                                                ! v
    CALL Initialise_model_vars (Model_data_array1(2), .TRUE.)                                                ! v

    DO loop = 1, 6
      PRINT *, 'Interpolation in 3d (2 time steps, quantity', loop, ')'
      WRITE (12,*) '--- Interpolation in 3d (2 time steps, quantity', loop, ') ---'

      ! Act with forward operator
      ObValue  = Interpolate3D (2, Model_data_array1(1:2),          &
                                dims % longs_u(4:5),                &
                                dims % half_levs(8:9),              &
                                AxisValuest(1:2),                   &
                                4, 8,                               &
                                dims % longs_u(4) + 0.004,          &
                                dims % half_levs(8) + 45.0,         &
                                AxisValuest(1) + 1.4,               &
                                loop )                                     ! Av
      ! Compute the LHS
      LHS      = ObValue * ObValue                                         ! (A v)^T A v

      ! Act with adjoint operator
      ! Initialize model data structure
      CALL Initialise_model_vars (Model_data_array2(1), .FALSE.)
      CALL Initialise_model_vars (Model_data_array2(2), .FALSE.)
      CALL Interpolate3D_adj (ObValue, 2, Model_data_array2(1:2),   &
                              dims % longs_u(4:5),                  &
                              dims % half_levs(8:9),                &
                              AxisValuest(1:2),                     &
                              4, 8,                                 &
                              dims % longs_u(4) + 0.004,            &
                              dims % half_levs(8) + 45.0,           &
                              AxisValuest(1) + 1.4,                 &
                              loop )                                       ! A^T A v
      ! Compute the RHS
      RHS      = InnerProdModelSpace (Model_data_array1(1),         &
                                      Model_data_array2(1),         &
                                      ignore_halos=.TRUE.) +        &
                 InnerProdModelSpace (Model_data_array1(2),         &
                                      Model_data_array2(2),         &
                                      ignore_halos=.TRUE.)                 ! (A^T A v)^T v
      WRITE (12,*) 'LHS = ', LHS
      WRITE (12,*) 'RHS = ', RHS
    END DO


    ! --------------------------------------------------------------------------
    ! A whole observations file (no time evolution in perturbations)
    PRINT *, 'A whole obs file (no time evolution in perts)'
    WRITE (12,*) '--- A whole obs file (no time evolution in perts) ---'

    ! The observations file
    Obs_filename = TRIM(datadir_Obs) // '/' // TRIM(Obs_file)
    PRINT *, 'Observations file ', TRIM(Obs_filename)

    ! Initialise the observation structure (a linked list)
    NULLIFY(Observations)

    ! Read-in the observations (for thorough test, make sure that file has observations
    ! of all observation types)
    CALL Read_Obs ( Observations, Obs_filename, dt, timestepsreq, dt_da, DAtimestepsreq, maxtime )

    !Allocate the array of data of DA time steps
    ALLOCATE (LSfc(0:DAtimestepsreq))

    ! The initial conditions of the linearisation stte
    LS_filename = TRIM(datadirABCfcs) // '/' // TRIM(LS_file)
    PRINT *, 'The file containing the LS: ', TRIM(LS_filename)

    ! Set state to zero
    CALL Initialise_model_vars (LSfc(0), .FALSE.)
    CALL Initialise_dims (dims)

    ! Read in initial conditions
    PRINT*, 'Reading in initial conditions ...'
    CALL Read_state_2d (LS_filename, LSfc(0), dims, 1)
    PRINT*, '-- done'

    ! Set some commonly-used constants
    CALL Set_ht_dep_cons (dims)
    PRINT *, 'About to run model for ', timestepsreq, ' time steps'
    CALL ABC_NL_ModelDriver_DA (LSfc(0:DAtimestepsreq), dims, timestepsreq, DAtimestepsreq)
    PRINT *, ' -- done'

    ! Run the non-linear observation operator
    PRINT *, 'Running the non-linear obs operator'
    CALL ModelObservations ( DAtimestepsreq, LSfc(0:DAtimestepsreq), &
                             dims, dt_da, 0, Observations, .FALSE. )
    PRINT *, ' -- done'

    ! Perturbations
    !Allocate the array of data of DA time steps
    ALLOCATE (ABC1(0:DAtimestepsreq))

    Pert_filename = TRIM(datadirABCperts) // '/' // TRIM(Pert_file)
    PRINT *, 'The file containing the perts: ', TRIM(Pert_filename)
    PRINT*, 'Reading in perturbations ...'
    CALL Read_state_2d (Pert_filename, ABC1(0), dims, 1)
    PRINT*, '-- done'

    ! Copy pert to all DA times
    DO loop = 1, DAtimestepsreq
      ABC1(loop) = ABC1(0)
    END DO

    ! Compute the model observations (linear operator)
    CALL ModelObservations_linear ( DAtimestepsreq, LSfc(0:DAtimestepsreq), &
                                    ABC1(0:DAtimestepsreq),                 &
                                    dims, dt_da, 0, Observations )         ! Av
    ! Compute the LHS
    LHS = InnerProduct_obsself (Observations, .FALSE., 1)                  ! (A v)^T A v

    ! Act with adjoint operator
    ! Allocate the array of data of DA time steps
    ALLOCATE (ABC2(0:DAtimestepsreq))
    DO loop = 0, DAtimestepsreq
      CALL Initialise_model_vars (ABC2(loop), .FALSE.)
    END DO

    CALL ModelObservations_adj ( DAtimestepsreq, LSfc(0:DAtimestepsreq), &
                                 ABC2(0:DAtimestepsreq),                 &
                                 dims, dt_da, 0, Observations, 2 )         ! A^T Av
    ! Compute the RHS
    RHS = 0.0
    DO loop = 0, DAtimestepsreq
      RHS = RHS + InnerProdModelSpace (ABC2(loop),            &
                                       ABC1(loop),            &
                                       ignore_halos=.FALSE.)               ! (A^T A v)^T v
    END DO

    WRITE (12,*) 'LHS = ', LHS
    WRITE (12,*) 'RHS = ', RHS

    ! Output the observations (including model obs, etc.)
    ! The observations file
    Obs_filename = TRIM(datadirTestDA) // '/Obs_processed.dat'
    PRINT *, 'Observations file ', TRIM(Obs_filename)
    CALL Write_Obs (Obs_filename, maxtime, dt, timestepsreq, dt_da, DAtimestepsreq, &
                    Observations)


    ! Dealocate
    DEALLOCATE (LSfc, ABC1, ABC2)
    ! Deallocate the observation structure
    CALL DeAllocate_Obs (Observations)

  END IF


! =============================================================================

  IF (RunInvTests) THEN
    PRINT*, '===== Running inverse tests ====='
    WRITE (12,*) '===== Running inverse tests ====='

    ! Take Pert as a sample increment, and do U U^-1
    ! Read in Pert
    PRINT*, 'Reading in perturbation ...'
    CALL Read_state_2d (Pert_filename, Pert1, dims, 1)
    PRINT*, '-- done'

    ! --------------------------------------------------------------------------
    WRITE (12,*) '--- Matrix inversion ---'
    ! Use the bal_r - bal_r covariance matrix inside the CVT structure as a test
    testmatrix(1:nlevs,1:nlevs) = CVT % Cov_rbalrbal(1:nlevs,1:nlevs)
    CALL InverseSymMat (nlevs, testmatrix(1:nlevs,1:nlevs))

    ! This matrix should return as an approximation of its inverse.  Check this
    WRITE (12,*) '--------------------------------'
    WRITE (12,*) 'Here is the result of mat^-1 mat'
    check(1:nlevs,1:nlevs) = MATMUL(testmatrix(1:nlevs,1:nlevs), CVT % Cov_rbalrbal(1:nlevs,1:nlevs))
    DO z = 1, nlevs
      WRITE (12,*) check(z,1:nlevs)
    END DO
    WRITE (12,*) '--------------------------------'
    WRITE (12,*) 'Here is the result of mat mat^-1'
        check(1:nlevs,1:nlevs) = MATMUL(CVT % Cov_rbalrbal(1:nlevs,1:nlevs), testmatrix(1:nlevs,1:nlevs))
    DO z = 1, nlevs
      WRITE (12,*) check(z,1:nlevs)
    END DO
    WRITE (12,*) '--------------------------------'


    ! --------------------------------------------------------------------------
    WRITE (12,*) '--- Helmholtz ---'

    ! Act with inverse operator (compute psi and chi from u and v)
    CALL Helmholtz_inv (CV_data1 % v1(1:nlongs,1:nlevs),     &
                        CV_data1 % v2(1:nlongs,1:nlevs),     &
                        Pert1 % u(0:nlongs+1,1:nlevs),       &
                        Pert1 % v(0:nlongs+1,1:nlevs))
    CALL Boundaries_CV (CV_data1, set_1=.TRUE., set_2=.TRUE.)

    InvTestOutput_filename = TRIM(datadirTestDA) // '/uv2psi.nc'
    CALL Write_one_field (InvTestOutput_filename, nlongs, nlevs, &
                          CV_data1 % v1(1:nlongs,1:nlevs), 'uv2psi')
    InvTestOutput_filename = TRIM(datadirTestDA) // '/uv2chi.nc'
    CALL Write_one_field (InvTestOutput_filename, nlongs, nlevs, &
                          CV_data1 % v2(1:nlongs,1:nlevs), 'uv2chi')

    ! Act with forward operator (compute u and v from psi and chi)
    CALL Helmholtz (CV_data1 % v1(0:nlongs+1,1:nlevs),     &
                    CV_data1 % v2(0:nlongs+1,1:nlevs),     &
                    Pert2 % u(1:nlongs,1:nlevs),           &
                    Pert2 % v(1:nlongs,1:nlevs))
    CALL Boundaries (Pert2, set_u=.TRUE., set_v=.TRUE.)

    InvTestOutput_filename = TRIM(datadirTestDA) // '/psichi2u.nc'
    CALL Write_one_field (InvTestOutput_filename, nlongs, nlevs, &
                          Pert2 % u(1:nlongs,1:nlevs), 'psichi2u')
    InvTestOutput_filename = TRIM(datadirTestDA) // '/psichi2v.nc'
    CALL Write_one_field (InvTestOutput_filename, nlongs, nlevs, &
                          Pert2 % v(1:nlongs,1:nlevs), 'psichi2v')


    ! --------------------------------------------------------------------------
    WRITE (12,*) '--- Fourier transform ---'

    ! Act with FFT operator
    CALL fft_real2spec (Pert1 % r(1:nlongs,1:nlevs),     &
                        CV_data1 % v3(1:nlongs,1:nlevs))
    CALL Boundaries_CV (CV_data1, set_3=.TRUE.)

    ! Act with the inverse FFT operator
    CALL fft_spec2real (CV_data1 % v3(1:nlongs,1:nlevs), &
                        Pert2 % r(1:nlongs,1:nlevs))

    InvTestOutput_filename = TRIM(datadirTestDA) // '/fft_test.nc'
    CALL Write_one_field (InvTestOutput_filename, nlongs, nlevs, &
                          Pert2 % r(1:nlongs,1:nlevs), 'r')



    ! --------------------------------------------------------------------------
    WRITE (12,*) '--- U_p ---'

    ! Act with inverse operator
    CALL U_p_inv (LS, CV_data1, Pert1, CVT % CVT_order,                                   &
                  CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab, &
                  CVT % CVT_param_opt_reg,                                                &
                  CVT % Regression(1:nlevs,1:nlevs), dims,                                &
                  .TRUE.,                                                                 &
                  datadirTestDA)

    ! Act with the forward operator
    CALL U_p (LS, CV_data1, Pert2, CVT % CVT_order,                                       &
              CVT % CVT_param_opt_gb, CVT % CVT_param_opt_hb, CVT % CVT_param_opt_ab,     &
              CVT % CVT_param_opt_reg,                                                    &
              CVT % Regression(1:nlevs,1:nlevs), dims)

    ! Output the result to check identical to LS
    InvTestOutput_filename = TRIM(datadirTestDA) // '/UpInvTest.nc'
    CALL Write_state_2d (InvTestOutput_filename, Pert2, dims, 1, 0, 0)


  END IF


  CLOSE (12)


END PROGRAM Master_TestSuite
