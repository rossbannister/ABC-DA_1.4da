#!/bin/bash
#set -x
# ========================================================
# To run an ABCvn1.4da to generate a forecast ensemble by running
# data assimilation (random background and random observations)
# followed by a forecast
# Author, Ross Bannister
# 05/08/2018
# ========================================================
# This script assumes that the following already exist:
#   CVT file, truth file, observation network specification file
# ========================================================


# ===== USER-SPECIFIED VARIABLES =========================

# Do a dummy run? (1 (yes) or 0 (no))
DUMMY=0

# Delete work files after script has run? (1 (yes) or 0 (no))
DELETE_WORK_FILES=1

# The directory at the base of this cycling (output data - MUST BE FULL PATH)
BASE_DIR=/home/basedir

# The starting member number
FIRST_MEMBER=401

# The finishing member number
LAST_MEMBER=512

# The number of data assimilation cycles
N_CYCLES=6

# The number of outer loops
N_OUTER_LOOPS=1

# The maximum number of inner loops
N_INNER_LOOPS_MAX=100

# Data assimilation window lenth (seconds)
DA_WINDOW=3600.0

# The directory containing the file that specifies which observations to make (input data only)
OBS_NETWORK_DIR=$BASE_DIR/../Master_MakeBgObs_ObsNetwork

# The directory containing the true state at the initial time (input data only)
INITIAL_TRUTH_DIR=$BASE_DIR/../da_cycle_0030/Obs+Truth

# The directory containing the CVT file (describing the background error covariances) (input data only)
CVT_DIR=$BASE_DIR/../../MasterCalibration

# The directory containing the main source files
CODE_DIR=/home/ABC_vn1.4da/src

# Set the Var assim type 3=3DVar, 35=3D-FGAT, 4=4DVar
VAR_TYPE=35


# Set the model parameters
A=0.02
B=0.01
C=1.0E4


# ===== END OF USER-SPECIFIED VARIABLES ==================



M1=-1

echo "=============================================" > $BASE_DIR/SuiteOut
echo "Running MakeEns script" >> $BASE_DIR/SuiteOut
date >> $BASE_DIR/SuiteOut
echo "=============================================" >> $BASE_DIR/SuiteOut


MEMBERDIR_PREV=0
MEMBERDIR_PREV_PREV=0

# =====================================================================================
# LOOP ROUND THE REQUIRED ENSEMBLE MEMBERS
# =====================================================================================
for MEMBER in $(seq $FIRST_MEMBER $LAST_MEMBER)
do
  MEMBER_FORM=`printf "%04u" $MEMBER`
  echo "Member " $MEMBER_FORM
  echo "========================================" >> $BASE_DIR/SuiteOut
  echo "===== Proceeding to generate ensemble member" $MEMBER_FORM "=====" >> $BASE_DIR/SuiteOut


  # Set the name of the (N_OUTER_LOOPS+1)th LS state produced by the DA
  # (the last time state in this file is the background for the next cycle)
  let N_OUTER_LOOPSP1=$N_OUTER_LOOPS+1
  OUTER_LOOP_FORM=`printf "%03u" $N_OUTER_LOOPSP1`
  USUAL_BACKGROUND_FILE=LS_Oloop${OUTER_LOOP_FORM}_Iloop000.nc

  MEMBER_DIR=$BASE_DIR/Member$MEMBER_FORM

  # Make the directory for this member
  mkdir -p $MEMBER_DIR

  # Set if the background state at the start has already been computed from the truth (1 or 0)
  BACKGROUND_ALREADY_COMPUTED=0


  # =====================================================================================
  # START THE CYCLING OF DATA ASSIMILATIONS
  # =====================================================================================
  for CYCLE in $(seq 1 $N_CYCLES)
  do
    CYCLE_FORM=`printf "%04u" $CYCLE`
    echo "===== DATA ASSIMILATION CYCLE NUMBER" $CYCLE_FORM "=====" >> $BASE_DIR/SuiteOut
    CYCLE_DIR=$MEMBER_DIR/da_cycle_$CYCLE_FORM

    # Make the directory for this da cycle
    mkdir -p $CYCLE_DIR

    if [ $CYCLE == 1 ]; then
      # Create a background state by perturbing the truth
      # -------------------------------------------------
      echo "  First cycle - need a background state" >> $BASE_DIR/SuiteOut
      if [ $BACKGROUND_ALREADY_COMPUTED == 0 ]; then
        mkdir -p $CYCLE_DIR/InitBg
        # Generate a random seed for background generation
        RANDOM_SEED=`date +%N`
        echo "Random number seed for initial background creation" $RANDOM_SEED >> $BASE_DIR/SuiteOut

        # Generate the namelist (see Sect. 4.8 of documentation, under Generate_mode=3)
        cat > $CYCLE_DIR/InitBg/UserOptions.nl << EOF
&UserOptions
! Generate background
! ------------------------------
  Generate_mode = 3
  datadirABC_in = '$INITIAL_TRUTH_DIR'
  init_ABC_file = 'Truth.nc'
  datadirCVT    = '$CVT_DIR'
  CVT_file      = 'CVT.nc'
  datadir_Bg    = '.'
  Pert_file     = 'Bgerr.nc'
  Bg_file       = 'Bg.nc'
  random_seed   = $RANDOM_SEED
/
EOF

        # Go into this directory and run the code to generate the background
        cd $CYCLE_DIR/InitBg
        if [ $DUMMY == 0 ]; then
          echo "  Generating the first background state ..." >> $BASE_DIR/SuiteOut
          $CODE_DIR/Master_MakeBgObs.out > stdout 2> stderr
          echo "  ... done" >> $BASE_DIR/SuiteOut
        fi
      fi

      # Set where the truth is for generating the obs
      TRUTH_DIR=$INITIAL_TRUTH_DIR
      TRUTH_FILE=Truth.nc

      # Make sure that the first DA cycle can find the background file
      CYCLE_DIR_PREV=$CYCLE_DIR/InitBg
      BACK_FILE=Bg.nc

    else
      # We have moved beyond the first cycle
      echo "  Background state to come from previous cycle" >> $BASE_DIR/SuiteOut
      # Set where the truth is for generating the obs
      TRUTH_DIR=$CYCLE_DIR_PREV/Obs+Truth
      TRUTH_FILE=Truth.nc
      # Set the name of the background file
      BACK_FILE=$USUAL_BACKGROUND_FILE
    fi

    echo "  This cycle's dir :" $CYCLE_DIR >> $BASE_DIR/SuiteOut
    echo "  Truth            :" $TRUTH_DIR/$TRUTH_FILE >> $BASE_DIR/SuiteOut
    echo "  Background       :" $CYCLE_DIR_PREV/$BACK_FILE >> $BASE_DIR/SuiteOut


    # Make the observations for this cycle
    # ------------------------------------
    # This also generates the truth state over the assimilation window
    mkdir -p $CYCLE_DIR/Obs+Truth
    # Generate a random seed for background generation
    RANDOM_SEED=`date +%N`
    echo "Random number seed for observation generation" $RANDOM_SEED >> $BASE_DIR/SuiteOut

    cat > $CYCLE_DIR/Obs+Truth/UserOptions.nl << EOF
&UserOptions
! Generate observations
! ------------------------------
  Generate_mode   = 2
  datadir_ObsSpec = '$OBS_NETWORK_DIR'
  ObsSpec_file    = 'ObsSpec.dat'
  datadir_Obs     = '.'
  Obs_file        = 'Obs.dat'
  datadirABC_in   = '$TRUTH_DIR'
  init_ABC_file   = '$TRUTH_FILE'
  output_ABC_file = 'Truth.nc'
  dt_da           = 600.0
  t0              = 0
  runlength       = $DA_WINDOW
  A               = $A
  B               = $B
  C               = $C
  random_seed     = $RANDOM_SEED
/
EOF

    # Go into this directory and run the code to generate the observations
    cd $CYCLE_DIR/Obs+Truth
    if [ $DUMMY == 0 ]; then
      echo "  Generating the synthetic obs ..." >> $BASE_DIR/SuiteOut
      $CODE_DIR/Master_MakeBgObs.out > stdout 2> stderr
      echo "  ... done" >> $BASE_DIR/SuiteOut
    fi


    # Run the data assimilation for this cycle
    # ----------------------------------------
    # At the end this also produces the background for the next cycle
    cat > $CYCLE_DIR/UserOptions.nl << EOF
&UserOptions
! Run the data assimilation
! -------------------------
  Vartype          = $VAR_TYPE,          !3=3DVar, 35=3D-FGAT, 4=4DVar
  Hybrid_opt       = 1,                  !Type of hybrid (or if pure Var)
  datadir_Bg       = '$CYCLE_DIR_PREV',
  Bg_file          = '$BACK_FILE',
  datadirCVT       = '$CVT_DIR',
  CVT_file         = 'CVT.nc',
  datadir_Obs      = 'Obs+Truth',
  Obs_file         = 'Obs.dat',
  t0               = 0,                  !Time of start of this DA cycle
  N_outerloops     = $N_OUTER_LOOPS,
  N_innerloops_max = $N_INNER_LOOPS_MAX,
  crit_inner       = 0.0001,
  datadirAnal      = '.',
  anal_file        = 'Anal.nc',
  analinc_file     = 'AnalInc.nc',
  diagnostics_file = 'diagnostics.dat'
/
EOF

    # Go into this directory and run the code to do the data assimilation
    cd $CYCLE_DIR
    if [ $DUMMY == 0 ]; then
      echo "  Assimilation ..." >> $BASE_DIR/SuiteOut
      $CODE_DIR/Master_Assimilate.out > stdout 2> stderr
      echo "  ... done" >> $BASE_DIR/SuiteOut
    fi
    # Prepare for the next cycle
    CYCLE_DIR_PREV=$CYCLE_DIR

  done



  # Extract the last time from the forecast file of the last cycle to use as the ensemble member
  # This involves running the model with 0.0 runlength and 0 dumps
  cd ../
  cat > UserOptions.nl << EOF
&UserOptions
! Running ABC model
! ------------------------------
  datadirABC_in            = 'da_cycle_$CYCLE_FORM'
  init_ABC_file            = '$USUAL_BACKGROUND_FILE'
  datadirABC_out           = '..'
  output_ABC_file          = 'Member$MEMBER_FORM.nc'
  diagnostics_file         = 'ABC_Diagnostics.dat'
  runlength                = 0.0
  ndumps                   = 0
  dt                       = 1.0
  Lengthscale_diagnostics  = .FALSE.
  A                        = 0.02
  B                        = 0.01
  C                        = 1.0E4
  Adv_tracer               = .TRUE.
/
EOF

  if [ $DUMMY == 0 ]; then
    echo "  Extracting forecast state ..." >> $BASE_DIR/SuiteOut
    $CODE_DIR/Master_RunNLModel.out > stdout 2> stderr
    echo "  ... done" >> $BASE_DIR/SuiteOut
  fi

  # Delete the work files
  if [ $DELETE_WORK_FILES == 1 ]; then
    echo "Deleting work files" >> $BASE_DIR/SuiteOut
    rm -r $BASE_DIR/Member$MEMBER_FORM
  fi

  echo "Produced ensemble member:" $MEMBER_FORM "in file" Member$MEMBER_FORM.nc >> $BASE_DIR/SuiteOut
done



echo "Script finished" >> $BASE_DIR/SuiteOut
echo "Script finished"
