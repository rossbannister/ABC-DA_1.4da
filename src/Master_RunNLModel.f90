PROGRAM Master_RunNLModel

!*****************************************************
!*   Code to run the non-linear ABC model            *
!*                                                   *
!*   R. Petrie,    2.0:  10-06-2011                  *
!*   R. Petrie,    3.0:  30-07-2013                  *
!*   R. Bannister, 3.1:  30-07-2013                  *
!*   R. Bannister, 1.4da 22-10-2017                  *
!*                                                   *
!*****************************************************


! Use Statements
!===============

USE DefConsTypes, ONLY :         &
    dims_type,                   &
    ABC_type,                    &
    datadirABC_in,               &
    datadirABC_out,              &
    init_ABC_file,               &
    output_ABC_file,             &
    ntimesteps,                  &
    ndumps,                      &
    diagnostics_file


IMPLICIT NONE

! Declare variables
!==========================
TYPE(dims_type)          :: dims
TYPE(ABC_type)           :: ABC_data

CHARACTER(LEN=320)       :: ABC_init_file, ABC_output_file, ABC_diags_file




PRINT*, '*************************************************************************'
PRINT*, 'Running Master_RunNLModel'
PRINT*, '*************************************************************************'

  ! Read namelist
  CALL SetOptions

  ABC_init_file   = TRIM(datadirABC_in)  // '/' // TRIM(init_ABC_file)
  ABC_output_file = TRIM(datadirABC_out) // '/' // TRIM(output_ABC_file)
  ABC_diags_file  = TRIM(datadirABC_out) // '/' // TRIM(diagnostics_file)

  ! Set state to zero
  CALL Initialise_model_vars (ABC_data, .FALSE.)
  CALL Initialise_dims (dims)

  ! Read in preprocessed UM data store in ABC_data
  PRINT*, 'Reading in processed data ...'
  CALL Read_state_2d (ABC_init_file, ABC_data, dims, -1)
  PRINT*, '-- done'

  ! Set some commonly-used constants
  CALL Set_ht_dep_cons (dims)

  CALL ABC_NL_ModelDriver ( ABC_data, dims, ntimesteps, ndumps,      &
                            ABC_output_file, ABC_diags_file )

END PROGRAM Master_RunNLModel
