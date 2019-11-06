#================================================================================
# Compile and link options for ubuntu
#===================================================================================
CompileOpts=-O3
#CompileOpts=
LinkOpts=-L/usr/local/lib -lnetcdff -lfftpack5.1d
LinkOptsEig=-L/usr/local/lib -lnetcdff -lfftpack5.1d -llapack -ltmglib -lrefblas

#================================================================================
# Compile and link options for Met Dept cluster
#===================================================================================
#CompileOpts=-O3
#CompileOpts=
#LinkOpts=-L/opt/software/modules/Met -L/home/ross/Lib/lib -lnetcdff -lfftpack5.1d



# Linking stage
Master_PrepareABC_InitState.out: Master_PrepareABC_InitState.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Write_state_2d.o BoundaryMod.o Boundaries.o Functions.o Energy.o Calc_geost.o Calc_hydro.o Initialise_CVT.o Initialise_ObsSpec.o Effective_buoyancy.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o
	f95 -o Master_PrepareABC_InitState.out Master_PrepareABC_InitState.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Write_state_2d.o BoundaryMod.o Boundaries.o Functions.o Energy.o Calc_geost.o Calc_hydro.o Initialise_CVT.o Initialise_ObsSpec.o Effective_buoyancy.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o $(LinkOpts)

Master_RunNLModel.out: Master_RunNLModel.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o
	f95 -o Master_RunNLModel.out Master_RunNLModel.o DefConsTypes.o SetOptions.o Initialise_dims.o Initialise_model_vars.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o ABC_NL_model.o Initialise_Averages.o Lscales_from_fft.o Magnitude_rms.o Boundaries.o Calc_geost.o Calc_hydro.o Energy.o Write_state_2d.o Functions.o Initialise_CVT.o Initialise_ObsSpec.o $(LinkOpts)

Master_Linear_Analysis.out: Master_Linear_Analysis.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o
	f95 -o Master_Linear_Analysis.out Master_Linear_Analysis.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o $(LinkOptsEig)

Master_TestSuite.out: Master_TestSuite.o DefConsTypes.o SetOptions.o Read_state_2d.o Set_ht_dep_cons.o Initialise_model_vars.o Functions.o Boundaries.o Boundaries_adj.o Initialise_CVs.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_adj.o U_h.o U_h_adj.o U_trans.o U_trans_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Write_state_2d.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Write_one_field.o Initialise_CVT.o Initialise_ObsSpec.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o Initialise_dims.o ABC_NL_ModelDriver_DA.o ABC_NL_model.o Initialise_Averages.o ModelObservations.o Write_Obs.o DeAllocate_Obs.o
	f95 -o Master_TestSuite.out Master_TestSuite.o DefConsTypes.o SetOptions.o Read_state_2d.o Set_ht_dep_cons.o Initialise_model_vars.o Functions.o Boundaries.o Boundaries_adj.o Initialise_CVs.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_adj.o U_h.o U_h_adj.o U_trans.o U_trans_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Write_state_2d.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Write_one_field.o Initialise_CVT.o Initialise_ObsSpec.o  Read_Obs.o ModelObservations_linear.o ModelObservations_adj.o Initialise_dims.o ABC_NL_ModelDriver_DA.o ABC_NL_model.o Initialise_Averages.o ModelObservations.o Write_Obs.o DeAllocate_Obs.o $(LinkOptsEig)

Master_Calibration.out: Master_Calibration.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Boundaries.o Write_state_2d.o ABC_NL_model.o Functions.o BoundaryMod.o Initialise_Averages.o Read_state_2d.o Helmholtz_inv.o LinearBal_r.o Calc_vert_cov1.o InverseSymMat.o Boundaries_CV.o Add_model_vars.o Div_model_cons.o Write_CV.o Subtract_model_vars.o U_p_inv.o Anbalw.o HydroBal_b.o Write_one_field.o Read_Covs.o Add_CVs.o Div_CV_cons.o Subtract_CVs.o CVT_Calibration_stddevs.o CVT_Calibration_vertcovs.o VertEigens.o CVT_Calibration_horizcovs.o Write_Covs.o Read_CV.o U_stddev_inv.o U_h_inv.o U_v_inv.o Initialise_CVT.o Initialise_CVs.o  fft_real2spec.o Initialise_ObsSpec.o Smooth.o Ensure_correlation_matrix.o
	f95 -o Master_Calibration.out Master_Calibration.o DefConsTypes.o SetOptions.o Initialise_um_data.o Initialise_dims.o Initialise_model_vars.o Read_um_data_2d.o Set_grid.o Set_ht_dep_cons.o Process_um_data.o Boundaries.o Write_state_2d.o ABC_NL_model.o Functions.o BoundaryMod.o Initialise_Averages.o Read_state_2d.o Helmholtz_inv.o LinearBal_r.o Calc_vert_cov1.o InverseSymMat.o Boundaries_CV.o Add_model_vars.o Div_model_cons.o Write_CV.o Subtract_model_vars.o U_p_inv.o Anbalw.o HydroBal_b.o Write_one_field.o Read_Covs.o Add_CVs.o Div_CV_cons.o Subtract_CVs.o CVT_Calibration_stddevs.o CVT_Calibration_vertcovs.o VertEigens.o CVT_Calibration_horizcovs.o Write_Covs.o Read_CV.o U_stddev_inv.o U_h_inv.o U_v_inv.o Initialise_CVT.o Initialise_CVs.o  fft_real2spec.o Initialise_ObsSpec.o Smooth.o Ensure_correlation_matrix.o $(LinkOptsEig)

Master_MakeBgObs.out: Master_MakeBgObs.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o DeAllocate_Obs.o Initialise_model_vars.o Initialise_dims.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver_DA.o ModelObservations.o Functions.o Boundaries.o ABC_NL_model.o Initialise_Averages.o Write_state_2d.o Read_Covs.o Initialise_CVs.o U_trans.o Add_model_vars.o U_h.o U_v.o U_p.o U_stddev.o fft_spec2real.o Helmholtz.o LinearBal_r.o HydroBal_b.o Anbalw.o Boundaries_CV.o Write_Obs.o Calc_hydro.o Calc_geost.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o
	f95 -o Master_MakeBgObs.out Master_MakeBgObs.o DefConsTypes.o SetOptions.o Initialise_CVT.o Initialise_ObsSpec.o DeAllocate_Obs.o Initialise_model_vars.o Initialise_dims.o Read_state_2d.o Set_ht_dep_cons.o ABC_NL_ModelDriver_DA.o ModelObservations.o Functions.o Boundaries.o ABC_NL_model.o Initialise_Averages.o Write_state_2d.o Read_Covs.o Initialise_CVs.o U_trans.o Add_model_vars.o U_h.o U_v.o U_p.o U_stddev.o fft_spec2real.o Helmholtz.o LinearBal_r.o HydroBal_b.o Anbalw.o Boundaries_CV.o Write_Obs.o Calc_hydro.o Calc_geost.o Calc_vert_mom_source.o Calc_horiz_div.o Calc_horiz_vort.o Effective_buoyancy.o $(LinkOpts)

Master_ImpliedCov.out: Master_ImpliedCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_trans_adj.o U_trans.o Write_state_2d.o U_v_adj.o U_h_adj.o U_p_adj.o U_stddev.o Boundaries_adj.o Anbalw_adj.o HydroBal_b_adj.o LinearBal_r_adj.o Helmholtz_adj.o U_v.o U_h.o U_p.o Boundaries.o Anbalw.o HydroBal_b.o LinearBal_r.o Helmholtz.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o
	f95 -o Master_ImpliedCov.out Master_ImpliedCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_trans_adj.o U_trans.o Write_state_2d.o U_v_adj.o U_h_adj.o U_p_adj.o U_stddev.o Boundaries_adj.o Anbalw_adj.o HydroBal_b_adj.o LinearBal_r_adj.o Helmholtz_adj.o U_v.o U_h.o U_p.o Boundaries.o Anbalw.o HydroBal_b.o LinearBal_r.o Helmholtz.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o $(LinkOpts)

Master_ImpliedTests.out: Master_ImpliedTests.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_v_adj.o U_h_adj.o Boundaries_adj.o U_v.o U_h.o Boundaries.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o
	f95 -o Master_ImpliedTests.out Master_ImpliedTests.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_CVs.o U_v_adj.o U_h_adj.o Boundaries_adj.o U_v.o U_h.o Boundaries.o fft_real2spec.o fft_spec2real.o Initialise_CVT.o Initialise_ObsSpec.o Functions.o Boundaries_CV.o Set_ht_dep_cons.o Read_Covs.o Initialise_model_vars.o Write_CV.o Boundaries_CV_adj.o $(LinkOpts)

Master_RawCov.out: Master_RawCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Initialise_CVT.o Initialise_ObsSpec.o Boundaries.o
	f95 -o Master_RawCov.out Master_RawCov.o DefConsTypes.o SetOptions.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Initialise_CVT.o Initialise_ObsSpec.o Boundaries.o $(LinkOpts)

Master_Assimilate.out: Master_Assimilate.o DefConsTypes.o SetOptions.o PenAndGrad.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Functions.o Boundaries.o Boundaries_adj.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_adj.o U_h.o U_h_adj.o U_trans.o U_trans_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o DeAllocate_Obs.o Initialise_CVs.o Subtract_model_vars.o U_trans_inv.o Minus_CVs.o Add_pert_CVs.o Add_model_vars.o Initialise_ObsSpec.o ABC_NL_ModelDriver_DA.o ModelObservations_ZeroPert.o Add_CVs.o Subtract_CVs.o Write_one_field.o U_stddev_inv.o U_v_inv.o U_h_inv.o Initialise_CVT.o ModelObservations.o ABC_NL_model.o Initialise_Averages.o Write_CV.o Write_Obs.o Energy.o Calc_geost.o Calc_hydro.o
	f95 -o Master_Assimilate.out Master_Assimilate.o DefConsTypes.o SetOptions.o PenAndGrad.o Read_state_2d.o Initialise_model_vars.o Set_ht_dep_cons.o Write_state_2d.o Functions.o Boundaries.o Boundaries_adj.o Boundaries_CV.o Boundaries_CV_adj.o LinearBal_r.o LinearBal_r_adj.o Anbalw.o Anbalw_adj.o Helmholtz.o Helmholtz_adj.o HydroBal_b.o HydroBal_b_adj.o U_p.o U_p_adj.o U_v.o U_v_adj.o U_h.o U_h_adj.o U_trans.o U_trans_adj.o fft_real2spec.o fft_spec2real.o U_stddev.o U_p_inv.o Helmholtz_inv.o Read_Covs.o InverseSymMat.o Read_Obs.o ModelObservations_linear.o  ModelObservations_adj.o DeAllocate_Obs.o Initialise_CVs.o Subtract_model_vars.o U_trans_inv.o Minus_CVs.o Add_pert_CVs.o Add_model_vars.o Initialise_ObsSpec.o ABC_NL_ModelDriver_DA.o ModelObservations_ZeroPert.o Add_CVs.o Subtract_CVs.o Write_one_field.o U_stddev_inv.o U_v_inv.o U_h_inv.o Initialise_CVT.o ModelObservations.o ABC_NL_model.o Initialise_Averages.o Write_CV.o Write_Obs.o Energy.o Calc_geost.o Calc_hydro.o $(LinkOptsEig)


# Compiling stage
Master_PrepareABC_InitState.o: Master_PrepareABC_InitState.f90 DefConsTypes.o
	f95 -c Master_PrepareABC_InitState.f90 $(CompileOpts)

DefConsTypes.o: DefConsTypes.f90
	f95 -c DefConsTypes.f90 $(CompileOpts)

SetOptions.o: SetOptions.f90 DefConsTypes.o
	f95 -c SetOptions.f90 $(CompileOpts)

Initialise_um_data.o: Initialise_um_data.f90 DefConsTypes.o
	f95 -c Initialise_um_data.f90 $(CompileOpts)

Initialise_dims.o: Initialise_dims.f90 DefConsTypes.o
	f95 -c Initialise_dims.f90 $(CompileOpts)

Initialise_model_vars.o: Initialise_model_vars.f90 DefConsTypes.o
	f95 -c Initialise_model_vars.f90 $(CompileOpts)

Read_um_data_2d.o: Read_um_data_2d.f90 DefConsTypes.o
	f95 -c Read_um_data_2d.f90 $(CompileOpts)

Set_grid.o: Set_grid.f90 DefConsTypes.o
	f95 -c Set_grid.f90 $(CompileOpts)

Set_ht_dep_cons.o: Set_ht_dep_cons.f90 DefConsTypes.o
	f95 -c Set_ht_dep_cons.f90 $(CompileOpts)

Process_um_data.o: Process_um_data.f90 DefConsTypes.o
	f95 -c Process_um_data.f90 $(CompileOpts)

Write_state_2d.o: Write_state_2d.f90 DefConsTypes.o
	f95 -c Write_state_2d.f90 $(CompileOpts)

Read_state_2d.o: Read_state_2d.f90 DefConsTypes.o
	f95 -c Read_state_2d.f90 $(CompileOpts)

BoundaryMod.o: BoundaryMod.f90 DefConsTypes.o
	f95 -c BoundaryMod.f90 $(CompileOpts)

Boundaries.o: Boundaries.f90 DefConsTypes.o
	f95 -c Boundaries.f90 $(CompileOpts)

Boundaries_adj.o: Boundaries_adj.f90 DefConsTypes.o
	f95 -c Boundaries_adj.f90 $(CompileOpts)

Functions.o: Functions.f90 DefConsTypes.o
	f95 -c Functions.f90 $(CompileOpts)

Energy.o: Energy.f90 DefConsTypes.o
	f95 -c Energy.f90 $(CompileOpts)

Calc_geost.o: Calc_geost.f90 DefConsTypes.o
	f95 -c Calc_geost.f90 $(CompileOpts)

Calc_hydro.o: Calc_hydro.f90 DefConsTypes.o
	f95 -c Calc_hydro.f90 $(CompileOpts)

Master_RunNLModel.o: Master_RunNLModel.f90 DefConsTypes.o
	f95 -c Master_RunNLModel.f90 $(CompileOpts)

ABC_NL_ModelDriver.o: ABC_NL_ModelDriver.f90 DefConsTypes.o
	f95 -c ABC_NL_ModelDriver.f90 $(CompileOpts)

Calc_vert_mom_source.o: Calc_vert_mom_source.f90 DefConsTypes.o
	f95 -c Calc_vert_mom_source.f90 $(CompileOpts)

Calc_horiz_div.o: Calc_horiz_div.f90 DefConsTypes.o
	f95 -c Calc_horiz_div.f90 $(CompileOpts)

Calc_horiz_vort.o: Calc_horiz_vort.f90 DefConsTypes.o
	f95 -c Calc_horiz_vort.f90 $(CompileOpts)

Effective_buoyancy.o: Effective_buoyancy.f90 DefConsTypes.o
	f95 -c Effective_buoyancy.f90 $(CompileOpts)

ABC_NL_model.o: ABC_NL_model.f90 DefConsTypes.o
	f95 -c ABC_NL_model.f90 $(CompileOpts)

Initialise_Averages.o: Initialise_Averages.f90 DefConsTypes.o
	f95 -c Initialise_Averages.f90 $(CompileOpts)

Lscales_from_fft.o: Lscales_from_fft.f90 DefConsTypes.o
	f95 -c Lscales_from_fft.f90 $(CompileOpts)

Magnitude_rms.o: Magnitude_rms.f90 DefConsTypes.o
	f95 -c Magnitude_rms.f90 $(CompileOpts)

Master_Linear_Analysis.o: Master_Linear_Analysis.f90 DefConsTypes.o
	f95 -c Master_Linear_Analysis.f90 $(CompileOpts)

Master_TestSuite.o: Master_TestSuite.f90 DefConsTypes.o
	f95 -c Master_TestSuite.f90 $(CompileOpts)

LinearBal_r.o: LinearBal_r.f90 DefConsTypes.o
	f95 -c LinearBal_r.f90 $(CompileOpts)

LinearBal_r_adj.o: LinearBal_r_adj.f90 DefConsTypes.o
	f95 -c LinearBal_r_adj.f90 $(CompileOpts)

Boundaries_CV.o: Boundaries_CV.f90 DefConsTypes.o
	f95 -c Boundaries_CV.f90 $(CompileOpts)

Boundaries_CV_adj.o: Boundaries_CV_adj.f90 DefConsTypes.o
	f95 -c Boundaries_CV_adj.f90 $(CompileOpts)

Initialise_CVs.o: Initialise_CVs.f90 DefConsTypes.o
	f95 -c Initialise_CVs.f90 $(CompileOpts)

Anbalw.o: Anbalw.f90 DefConsTypes.o
	f95 -c Anbalw.f90 $(CompileOpts)

Anbalw_adj.o: Anbalw_adj.f90 DefConsTypes.o
	f95 -c Anbalw_adj.f90 $(CompileOpts)

Helmholtz.o: Helmholtz.f90 DefConsTypes.o
	f95 -c Helmholtz.f90 $(CompileOpts)

Helmholtz_adj.o: Helmholtz_adj.f90 DefConsTypes.o
	f95 -c Helmholtz_adj.f90 $(CompileOpts)

Helmholtz_inv.o: Helmholtz_inv.f90 DefConsTypes.o
	f95 -c Helmholtz_inv.f90 $(CompileOpts)

HydroBal_b.o: HydroBal_b.f90 DefConsTypes.o
	f95 -c HydroBal_b.f90 $(CompileOpts)

HydroBal_b_adj.o: HydroBal_b_adj.f90 DefConsTypes.o
	f95 -c HydroBal_b_adj.f90 $(CompileOpts)

U_p.o: U_p.f90 DefConsTypes.o
	f95 -c U_p.f90 $(CompileOpts)

U_p_adj.o: U_p_adj.f90 DefConsTypes.o
	f95 -c U_p_adj.f90 $(CompileOpts)

U_v.o: U_v.f90 DefConsTypes.o
	f95 -c U_v.f90 $(CompileOpts)

U_v_adj.o: U_v_adj.f90 DefConsTypes.o
	f95 -c U_v_adj.f90 $(CompileOpts)

U_h.o: U_h.f90 DefConsTypes.o
	f95 -c U_h.f90 $(CompileOpts)

U_h_adj.o: U_h_adj.f90 DefConsTypes.o
	f95 -c U_h_adj.f90 $(CompileOpts)

U_trans.o: U_trans.f90 DefConsTypes.o
	f95 -c U_trans.f90 $(CompileOpts)

U_trans_adj.o: U_trans_adj.f90 DefConsTypes.o
	f95 -c U_trans_adj.f90 $(CompileOpts)

fft_real2spec.o: fft_real2spec.f90 DefConsTypes.o
	f95 -c fft_real2spec.f90 $(CompileOpts)

fft_spec2real.o: fft_spec2real.f90 DefConsTypes.o
	f95 -c fft_spec2real.f90 $(CompileOpts)

U_stddev.o: U_stddev.f90 DefConsTypes.o
	f95 -c U_stddev.f90 $(CompileOpts)

U_p_inv.o: U_p_inv.f90 DefConsTypes.o
	f95 -c U_p_inv.f90 $(CompileOpts)

Master_Calibration.o: Master_Calibration.f90 DefConsTypes.o
	f95 -c Master_Calibration.f90 $(CompileOpts)

Calc_vert_cov1.o: Calc_vert_cov1.f90 DefConsTypes.o
	f95 -c Calc_vert_cov1.f90 $(CompileOpts)

InverseSymMat.o: InverseSymMat.f90 DefConsTypes.o
	f95 -c InverseSymMat.f90 $(CompileOpts)

Write_Covs.o: Write_Covs.f90 DefConsTypes.o
	f95 -c Write_Covs.f90 $(CompileOpts)

Add_model_vars.o: Add_model_vars.f90 DefConsTypes.o
	f95 -c Add_model_vars.f90 $(CompileOpts)

Div_model_cons.o: Div_model_cons.f90 DefConsTypes.o
	f95 -c Div_model_cons.f90 $(CompileOpts)

Write_CV.o: Write_CV.f90 DefConsTypes.o
	f95 -c Write_CV.f90 $(CompileOpts)

Subtract_model_vars.o: Subtract_model_vars.f90 DefConsTypes.o
	f95 -c Subtract_model_vars.f90 $(CompileOpts)

Write_one_field.o: Write_one_field.f90 DefConsTypes.o
	f95 -c Write_one_field.f90 $(CompileOpts)

Read_CV.o: Read_CV.f90 DefConsTypes.o
	f95 -c Read_CV.f90 $(CompileOpts)

Read_Covs.o: Read_Covs.f90 DefConsTypes.o
	f95 -c Read_Covs.f90 $(CompileOpts)

Add_CVs.o: Add_CVs.f90 DefConsTypes.o
	f95 -c Add_CVs.f90 $(CompileOpts)

Div_CV_cons.o: Div_CV_cons.f90 DefConsTypes.o
	f95 -c Div_CV_cons.f90 $(CompileOpts)

Subtract_CVs.o: Subtract_CVs.f90 DefConsTypes.o
	f95 -c Subtract_CVs.f90 $(CompileOpts)

CVT_Calibration_stddevs.o: CVT_Calibration_stddevs.f90 DefConsTypes.o
	f95 -c CVT_Calibration_stddevs.f90 $(CompileOpts)

CVT_Calibration_vertcovs.o: CVT_Calibration_vertcovs.f90 DefConsTypes.o
	f95 -c CVT_Calibration_vertcovs.f90 $(CompileOpts)

VertEigens.o: VertEigens.f90 DefConsTypes.o
	f95 -c VertEigens.f90 $(CompileOpts)

CVT_Calibration_horizcovs.o: CVT_Calibration_horizcovs.f90 DefConsTypes.o
	f95 -c CVT_Calibration_horizcovs.f90 $(CompileOpts)

Initialise_CVT.o: Initialise_CVT.f90 DefConsTypes.o
	f95 -c Initialise_CVT.f90 $(CompileOpts)

U_stddev_inv.o: U_stddev_inv.f90 DefConsTypes.o
	f95 -c U_stddev_inv.f90 $(CompileOpts)

U_h_inv.o: U_h_inv.f90 DefConsTypes.o
	f95 -c U_h_inv.f90 $(CompileOpts)

U_v_inv.o: U_v_inv.f90 DefConsTypes.o
	f95 -c U_v_inv.f90 $(CompileOpts)

Master_MakeBgObs.o: Master_MakeBgObs.f90 DefConsTypes.o
	f95 -c Master_MakeBgObs.f90 $(CompileOpts)

Initialise_ObsSpec.o: Initialise_ObsSpec.f90 DefConsTypes.o
	f95 -c Initialise_ObsSpec.f90 $(CompileOpts)

DeAllocate_Obs.o: DeAllocate_Obs.f90 DefConsTypes.o
	f95 -c DeAllocate_Obs.f90 $(CompileOpts)

ABC_NL_ModelDriver_DA.o: ABC_NL_ModelDriver_DA.f90 DefConsTypes.o
	f95 -c ABC_NL_ModelDriver_DA.f90 $(CompileOpts)

ModelObservations.o: ModelObservations.f90 DefConsTypes.o
	f95 -c ModelObservations.f90 $(CompileOpts)

ModelObservations_linear.o: ModelObservations_linear.f90 DefConsTypes.o
	f95 -c ModelObservations_linear.f90 $(CompileOpts)

ModelObservations_adj.o: ModelObservations_adj.f90 DefConsTypes.o
	f95 -c ModelObservations_adj.f90 $(CompileOpts)

Master_ImpliedCov.o: Master_ImpliedCov.f90 DefConsTypes.o
	f95 -c Master_ImpliedCov.f90 $(CompileOpts)

Master_ImpliedTests.o: Master_ImpliedTests.f90 DefConsTypes.o
	f95 -c Master_ImpliedTests.f90 $(CompileOpts)

Master_RawCov.o: Master_RawCov.f90 DefConsTypes.o
	f95 -c Master_RawCov.f90 $(CompileOpts)

Read_Obs.o: Read_Obs.f90 DefConsTypes.o
	f95 -c Read_Obs.f90 $(CompileOpts)

Write_Obs.o: Write_Obs.f90 DefConsTypes.o
	f95 -c Write_Obs.f90 $(CompileOpts)

Master_Assimilate.o: Master_Assimilate.f90 DefConsTypes.o
	f95 -c Master_Assimilate.f90 $(CompileOpts)

PenAndGrad.o: PenAndGrad.f90 DefConsTypes.o
	f95 -c PenAndGrad.f90 $(CompileOpts)

U_trans_inv.o: U_trans_inv.f90 DefConsTypes.o
	f95 -c U_trans_inv.f90 $(CompileOpts)

Minus_CVs.o: Minus_CVs.f90 DefConsTypes.o
	f95 -c Minus_CVs.f90 $(CompileOpts)

Add_pert_CVs.o: Add_pert_CVs.f90 DefConsTypes.o
	f95 -c Add_pert_CVs.f90 $(CompileOpts)

ModelObservations_ZeroPert.o: ModelObservations_ZeroPert.f90 DefConsTypes.o
	f95 -c ModelObservations_ZeroPert.f90 $(CompileOpts)

Smooth.o: Smooth.f90 DefConsTypes.o
	f95 -c Smooth.f90 $(CompileOpts)

Ensure_correlation_matrix.o: Ensure_correlation_matrix.f90 DefConsTypes.o
	f95 -c Ensure_correlation_matrix.f90 $(CompileOpts)

all: Master_PrepareABC_InitState.out Master_RunNLModel.out Master_Linear_Analysis.out Master_TestSuite.out Master_Calibration.out Master_MakeBgObs.out Master_ImpliedCov.out Master_RawCov.out Master_Assimilate.out

clean:
	rm -rf *.o *.out *.mod
