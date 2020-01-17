/* ==================================================================================
  Program to compute vertical lengthscales and Rossby radii as a function of vertical mode.
  For use with ABC-DA

  Language: C++

  Modification history
  --------------------
  12/09/19 New Code. Ross Bannister

  On Ubuntu machines
  g++ -I/usr/include RossbyRadius.cpp -o RossbyRadius.out -L/usr/lib -lnetcdf_c++ -lnetcdf -lfftw3 -lm

  Example usage
  ./RossbyRadius.out \
    /home/ross/DataAssim/RuthsModel/ABC_vn1.4da/Investigations/Orig_transform_order_smoothed_sigma/Exp+GB+HB-AB/MasterCalibration \
    /home/ross/DataAssim/RuthsModel/ABC_vn1.4da/examples/Master_Linear_Analysis


  =============================================================================== */

  #include <stdio.h>
  #include <math.h>
  #include <stdlib.h>
  #include <string.h>
  #include <netcdf.h>
  #include <fftw3.h>


/* -------------------------------------------------------------------------------
   Global constants
   ------------------------------------------------------------------------------- */

  const double  pi            = 3.14159265;            // pi
  const int     nlong         = 360;                   // number of longitudes
  const int     nlev          = 60;                    // number of levels
  const double  gridlength    = 1500.0;                // horizontal gridlength
  const double  g             = 9.81;                  // acceleration due to gravity

/* -------------------------------------------------------------------------------
   Global variables
   ------------------------------------------------------------------------------- */

/* -------------------------------------------------------------------------------
   Define the important data types
   ------------------------------------------------------------------------------- */



/* -------------------------------------------------------------------------------
   Function declarations
   ------------------------------------------------------------------------------- */

bool ProgramArguments (
       int   argument_count,           // in  Number of arguments entered
       char  **argument_list,          // in  The program arguments
       char  cvt_dir[256],             // out Name of directory containing CVT.nc
       char  lin_anal_dir[256]         // out Name of directory containing hori_grav_speed.dat
                      );

void Array_2d_double_create ( double ***Array,
                              int    xlen,
                              int    ylen );

void Array_2d_double_destroy ( double ***Array,
                               int    xlen);


bool read_vertmodes ( char    filename[],
                      double  *A,
                      double  *B,
                      double  *C,
                      double  *f,
                      double  *levels,
                      double  **VertModes,
                      int     nlev );

void EstimateVertLengths ( double **VertModes,
                           double VertLength[],
                           int    EffVertWn[],
                           double L,
                           int    nlev );

void EstimateRossbyRad ( double VertLength[],
                         double RossbyRad[],
                         int    EffRossbyWn[],
                         double f,
                         double g,
                         double L_horiz,
                         int    nlev );

void OutputResults (double **VertModes,
                    double VertLength[],
                    int    EffVertWn[],
                    double RossbyRad[],
                    int    EffRossbyWn[],
                    int    nlev,
                    char   filename_modes[],
                    char   filename_lengths[]);

void Read_horiz_gw_speeds ( double **Horiz_gw_speeds,
                            char   filename[],
                            int    nvertwns,
                            int    nhorizwns );

void Find_LRossby ( double **Horiz_gw_speeds,
                    int    nvertwns,
                    int    nhorizwns,
                    double L_horiz,
                    double f,
                    double RossbyRad[] );

/* -------------------------------------------------------------------------------
   Main part of the program
   ------------------------------------------------------------------------------- */
int main ( int   argument_count,
           char  **argument_list )


{ // Declare the main variables
  char                  cvt_dir[256], cvt_filename[280];
  char                  lin_anal_dir[256], lin_anal_filename[280];
  bool                  success, readok;
  double                A, B, C, f;
  double                **VertModes, levels[::nlev], VertLength[::nlev];
  int                   EffVertWn[::nlev];
  double                RossbyRad[::nlev];
  int                   EffRossbyWn[::nlev];
  double                **Horiz_gw_speeds;

  printf ("----------------------------------------------------------------------\n");
  printf ("Program to determine:\n  (a) vertical lengthscale\n  (b) Rossby radius\n");
  printf ("From vertical mode data\n");
  printf ("Ross Bannister, NCEO, 2019\n");


  // Retrieve the directory name containing the CVT data
  printf ("Extracting program arguments ...\n");
  success = ProgramArguments (argument_count,
                              argument_list,
                              cvt_dir,
                              lin_anal_dir);
  printf ("... done\n");

  if (success)
  { // Construct the total filenames
    sprintf (cvt_filename, "%s/CVT.nc", cvt_dir);
    sprintf (lin_anal_filename, "%s/hori_grav_speed.dat", lin_anal_dir);

    // Allocate space for the verical modes
    printf ("Allocating space for vertical mode data ...\n");
    Array_2d_double_create ( &VertModes,
                             ::nlev,
                             ::nlev );
    printf ("... done\n");

/*
    // Allocate space for the horizontal gravity wave speeds
    printf ("Allocating space for horizontal gravity wave data ...\n");
    Array_2d_double_create ( &Horiz_gw_speeds,
                             ::nlev,
                             ::nlong );
    printf ("... done\n");
*/

    // Read vertical modes, etc.
    printf ("Reading vertical modes from file %s ...\n", cvt_filename);
    readok = read_vertmodes (cvt_filename,
                             &A,
                             &B,
                             &C,
                             &f,
                             levels,
                             VertModes,
                             ::nlev);
    printf ("... done\n");

    if (readok)
    { printf ("A = %f\n", A);
      printf ("B = %f\n", B);
      printf ("C = %f\n", C);
      printf ("f = %f\n", f);
      printf ("L = %f\n", levels[::nlev-1]);

      // Estimate vertical lengthscales for each vertical mode
      printf ("Estimating vertical lengthscales by FFT ...\n");
      EstimateVertLengths ( VertModes,
                            VertLength,
                            EffVertWn,
                            levels[::nlev-1],
                            ::nlev );
      printf ("... done\n");

      // Estimate Rossby radius for each vertical mode
      printf ("Estimating Rossby radii ...\n");
      EstimateRossbyRad ( VertLength,
                          RossbyRad,
                          EffRossbyWn,
                          f,
                          ::g,
                          ::gridlength * ::nlong,
                          ::nlev );

/*
      // Read-in the horizontal gravity wave speeds
      printf ("Reading in the horizontal gravity wave speeds ...\n");
      Read_horiz_gw_speeds ( Horiz_gw_speeds,
                            lin_anal_filename,
                            ::nlev,
                            ::nlong );
      printf ("... done\n");

      // Find the most consistent Rossby radii
      printf ("Computing the most consistent Rossby radii ...\n");
      Find_LRossby ( Horiz_gw_speeds,
                     ::nlev,
                     ::nlong,
                     ::gridlength * ::nlong,
                     f,
                     RossbyRad );
*/

      // Output the computed data
      printf ("Outputting results ...\n");
      OutputResults ( VertModes,
                      VertLength,
                      EffVertWn,
                      RossbyRad,
                      EffRossbyWn,
                      ::nlev,
                      "VerticalModes.dat",
                      "VerticalModeCharacteristics.dat");
      printf ("... done\n");
    }

    // Deallocate
    Array_2d_double_destroy ( &VertModes,
                              ::nlev );

/*
    Array_2d_double_destroy ( &Horiz_gw_speeds,
                              ::nlev );
*/
  }
}




/* -------------------------------------------------------------------------------
   Function definitions
   ------------------------------------------------------------------------------- */

// -------------------------------------------------------------------------------/* -------------------------------------------------------------------------------
bool ProgramArguments (
       int   argument_count,           // in  Number of arguments entered
       char  **argument_list,          // in  The program arguments
       char  cvt_dir[256],             // out Name of directory of CVT file
       char  lin_anal_dir[256] )       // out Name of directory of hori_grav_speed.dat
// Extract key program arguments
{ // Declare local variables
  int    success, count;
  // Return 1 from subroutine if OK, 0 otherwise
  success      = true;
  count        = 0;

  count += 1;
  if (argument_count > count)
  { strcpy(cvt_dir, argument_list[count]);
    printf ("Name of cvt dir : %s\n", cvt_dir);
  }
  else
  { success = false;
    printf ("Argument %u is missing: please specify the name of directory of cvt file.\n", count);
  }

  count += 1;
  if (argument_count > count)
  { strcpy(lin_anal_dir, argument_list[count]);
    printf ("Name of linear analysis dir : %s\n", lin_anal_dir);
  }
  else
  { success = false;
    printf ("Argument %u is missing: please specify the name of directory of hori_grav_speed.dat file.\n", count);
  }

  return success;
}


// -------------------------------------------------------------------------------
void Array_2d_double_create ( double ***Array,
                              int    xlen,
                              int    ylen )
{ int j;
  *Array = new double *[xlen];
  for (j=0; j<xlen; j++)
  { (*Array)[j] = new double[ylen];
  }
}


// -------------------------------------------------------------------------------
void Array_2d_double_destroy ( double ***Array,
                               int    xlen)
{ int j;
  for (j=0; j<xlen; j++)
  { delete[] (*Array)[j];
  }
  delete[] (*Array);
}




// -------------------------------------------------------------------------------
bool read_vertmodes ( char    filename[],
                      double  *A,
                      double  *B,
                      double  *C,
                      double  *f,
                      double  *levels,
                      double  **VertModes,
                      int     nlev )

// Read-in A, B, C, f, and the vertical modes for unbalanced pressure
{ // Declare local variables
  bool   ok;
  int    ierr, ierr1, ierr2, ierr3, ierr4, ierr5, ncid;
  int    x, z;
  double lod[nlev];
  // netCDF-related variables
  int    varidlevel, varidA, varidB, varidC, varidf, varidvertmodes;
  size_t start1[1], count1[1];
  size_t start2[2], count2[2];


  ierr = nc_open (filename,
                  NC_NOWRITE,
                  &ncid);
  if (ierr != 0)
  { printf ("Error opening file %s\n%s\nExiting\n", filename, nc_strerror(ierr));
    exit(0);
  }

  // Get dimension variable ids
  ierr1 = nc_inq_varid (ncid,
                        "level",
                        &varidlevel);
  if (ierr1 != 0)
  { printf ("Error getting dimension variable ids\n%s\nExiting\n",
            nc_strerror(ierr1));
    exit(0);
  }

  // Get main field variable ids
  ierr1 = nc_inq_varid (ncid,
                        "A",
                        &varidA);
  ierr2 = nc_inq_varid (ncid,
                        "B",
                        &varidB);
  ierr3 = nc_inq_varid (ncid,
                        "C",
                        &varidC);
  ierr4 = nc_inq_varid (ncid,
                        "f",
                        &varidf);
  ierr5 = nc_inq_varid (ncid,
                        "vertmode3",
                        &varidvertmodes);

  if (ierr1 + ierr2 + ierr3 + ierr4 + ierr5 != 0)
  { printf ("Error getting main variable ids\n%s\n%s\n%s\n%s\n%s\nExiting\n",
            nc_strerror(ierr1),
            nc_strerror(ierr2),
            nc_strerror(ierr3),
            nc_strerror(ierr4),
            nc_strerror(ierr5));
    exit(0);
  }

  // Retrieve dimension variables
  // Levels
  start1[0] = 0;  count1[0] = nlev;
  ierr1     = nc_get_vara_double (ncid,
                                  varidlevel,
                                  start1,
                                  count1,
                                  levels);

  // Retrieve main variables
  // Retrieve parameters
  start1[0] = 0;  count1[0] = 1;
  ierr1     = nc_get_vara_double (ncid,
                                  varidA,
                                  start1,
                                  count1,
                                  lod);
  *A        = lod[0];
  ierr2     = nc_get_vara_double (ncid,
                                  varidB,
                                  start1,
                                  count1,
                                  lod);
  *B        = lod[0];
  ierr3     = nc_get_vara_double (ncid,
                                  varidC,
                                  start1,
                                  count1,
                                  lod);
  *C        = lod[0];
  ierr4     = nc_get_vara_double (ncid,
                                  varidf,
                                  start1,
                                  count1,
                                  lod);
  *f        = lod[0];
  if (ierr1 + ierr2 + ierr3 + ierr4 != 0)
  { printf ("Error getting parameter values\n%s\n%s\n%s\n%s\nExiting\n",
            nc_strerror(ierr1),
            nc_strerror(ierr2),
            nc_strerror(ierr3),
            nc_strerror(ierr4));
    exit(0);
  }

  // Retrieve vertical modes
  start2[1] = 0;     count2[1] = nlev;   // x
  start2[0] = 0;     count2[0] = 1;      // y
  for (z=0; z<nlev; z++)
  { start2[0] = z;
    ierr = nc_get_vara_double (ncid,
                               varidvertmodes,
                               start2,
                               count2,
                               lod);
    if (ierr != 0)
    { printf ("Error reading vertmodes\n%s\nExiting\n", nc_strerror(ierr));
      exit(0);
    }
    for (x=0; x<nlev; x++)
    { VertModes[x][z] = lod[x];
    }
  }


  //Close the netCDF file
  ierr = nc_close (ncid);
  if (ierr != 0)
  { printf ("Error closing output file\n%s\nExiting\n", nc_strerror(ierr));
    exit(0);
  }

  return true;
}


// -------------------------------------------------------------------------------
void EstimateVertLengths ( double **VertModes,
                           double VertLength[],
                           int    EffVertWn[],
                           double L,
                           int    nlev )
// Estimate the vertical lengthscales using info from the vertical modes
{ fftw_complex Spectral[nlev], Physical[nlev];
  fftw_plan    plan_ft;
  double       wn[nlev], avwn, norm, r, i, weight;
  int          k, wni, m, z;


  // Set-up 'plan' for forward and inverse FTs
  plan_ft  = fftw_plan_dft_1d (nlev, Physical, Spectral, FFTW_FORWARD, FFTW_ESTIMATE);

  // Set-up wavenumbers
  for (k=0; k<nlev; k++)
  { // Find wavenumber in the right zone
    wni   = (k < abs(k-nlev)) ? k : k-nlev;
    wn[k] = double(2*wni)*::pi/L;
  }

  // Set the imaginary part of the input structure to zero
  for (z=0; z<nlev; z++)
  { Physical[z][1] = 0.0;
  }

  // Loop over the vertical modes
  for (m=0; m<nlev; m++)
  { // Copy information from this vertical mode to the input structure of the fft
    for (z=0; z<nlev; z++)
    { Physical[z][0] = VertModes[z][m];
    }
    // Do Fourier transform of this
    fftw_execute(plan_ft);
    // Find the average wavenumber as though the FT is a PDF
    avwn = 0.0;
    norm = 0.0;
    for (k=0; k<nlev; k++)
    { r      = Spectral[k][0];
      i      = Spectral[k][1];
      weight = r*r + i*i;
      avwn  += fabs(wn[k]) * weight;
      norm  += weight;
    }
    avwn         /= norm;
    VertLength[m] = 2.0 * ::pi / avwn;
    EffVertWn[m]  = int(L * avwn / (2.0 * ::pi) + 0.5);
  }
}


// -------------------------------------------------------------------------------
void EstimateRossbyRad ( double VertLength[],
                         double RossbyRad[],
                         int    EffRossbyWn[],
                         double f,
                         double g,
                         double L_horiz,
                         int    nlev )
// Estimate the Rossby radii of each vertical mode
{ int m;

  for (m=0; m<nlev; m++)
  { RossbyRad[m]   = sqrt(g * VertLength[m]) / f;
    EffRossbyWn[m] = int(L_horiz / RossbyRad[m] + 0.5);
  }
}


// -------------------------------------------------------------------------------
void OutputResults (double **VertModes,
                    double VertLength[],
                    int    EffVertWn[],
                    double RossbyRad[],
                    int    EffRossbyWn[],
                    int    nlev,
                    char   filename_modes[],
                    char   filename_lengths[])
// Output mode pattern, lengthscale, and Rossby radius of each vert mode
{ FILE* file = NULL;
  int   m, z;

  // Output the mode patterns themselves
  file = fopen (filename_modes, "w");
  for (z=0; z<nlev; z++)
  { fprintf (file, "%02i  ", z);
    for (m=0; m<nlev; m++)
    { fprintf (file, "%f ", VertModes[z][m]);
    }
    fprintf (file, "\n");
  }
  fclose (file);

  // Output the vertical lengthscales, effective vert wavenumbers, and Rossby radii, and Rossby wn
  file = fopen (filename_lengths, "w");
  fprintf (file, "# index   Vert length   Eff vert wn   Rossby rad   Eff Rossby wn\n");
  for (m=0; m<nlev; m++)
  { fprintf (file, "%02i  %f  %i  %f  %i\n",
             m, VertLength[m], EffVertWn[m],
             RossbyRad[m], EffRossbyWn[m]);
  }
  fclose (file);
}


// -------------------------------------------------------------------------------
void Read_horiz_gw_speeds ( double **Horiz_gw_speeds,
                            char   filename[],
                            int    nvertwns,
                            int    nhorizwns )
// Read horizontal gravity wave speeds
{ FILE* file = NULL;
  int   m, k;
  float speed;

  printf ("Reading from file : %s\n", filename);
  file = fopen (filename, "r");
  for (m=0; m<nvertwns; m++)
  { for (k=0; k<nhorizwns; k++)
    { fscanf (file, "%f", &speed);
      Horiz_gw_speeds[m][k] = speed;
    }
    //printf ("%i : %f  %f  %f\n", m, Horiz_gw_speeds[m][0], Horiz_gw_speeds[m][1], Horiz_gw_speeds[m][nhorizwns-1]);
  }
  fclose (file);
}


// -------------------------------------------------------------------------------
void Find_LRossby ( double **Horiz_gw_speeds,
                    int    nvertwns,
                    int    nhorizwns,
                    double L_horiz,
                    double f,
                    double RossbyRad[] )
// Find the most consistent Rossby radius for each vertical wavenumber
{ int    m, k, bestk;
  double proposedRR[nhorizwns];
  double horizlengths[nhorizwns];
  double smallestdiff, diff;

  // Here is the strategy:
  //   - For each vertical wavenumber we have horiz gravity wave speeds as a fn of horiz wavenumber
  //   - Propose a Rossby radius for each gw speed
  //   - Compute a horiz lengthscale for each horiz wavenumber
  //   - Which ever proposed Rossby radius is closest to the horiz lengthscale IS taken as the actual Rossby radius

  // Compute the horizontal lengthscales (this is done once)
  for (k=1; k<nhorizwns; k++)
  { horizlengths[k] = L_horiz / float(k);
    // This is computed from the following equations
    // horizlen = 2pi / wavenumber
    // wavenumber = 2pi wn_index / L_horiz
    // note: wn_index = k
    // note: do not use k=0
  }

  // Loop over vertical wavenumber
  for (m=0; m<nvertwns; m++)
  { // Find the proposed Rossby radii for each
    for (k=1; k<nhorizwns; k++)
    { proposedRR[k] = Horiz_gw_speeds[m][k] / f;
    }
    // Which one is most consistent with the actual lengthscale?
    smallestdiff = 999999.0;
    bestk        = 0;
    for (k=1; k<nhorizwns; k++)
    { diff = fabs(proposedRR[k] - horizlengths[k]);
      if (diff < smallestdiff)
      { smallestdiff = diff;
        bestk        = k;
      }
    }
    if (k > 0)
    { // Found consistent value
      RossbyRad[m] = proposedRR[bestk];
    }
    else
    { // Failed to find value
      RossbyRad[m] = 0.0;
    }
    printf ("Vert mode %i : %f (diff %f)\n", m, RossbyRad[m], diff);
  }
}
