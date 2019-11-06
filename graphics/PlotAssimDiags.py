#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot assimilation diagnostics
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# The following are plotted
#   1. Background trajectory
#   2. Analysis trajectory
#   3. Analysis increment (t=0)
#   4. Propagation of the gradient wrt J back in time
#   5. Truth trajectory
#   6. Background error as a function of time
#   7. Analysis error as a function of time
#   8. Cost function with iteration
#   9. Energy with iteration, plus the true energy
#  10. Size of gradient of cost function with iteration
#  11. Imalance diagnostics with time (background, analysis, and truth)
#  12. Observations, background observations, analysis observations
#
# Ross Bannister, July 2018
# -------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os
import sys
from Routines4PlotAssimDiags import *

# Set variables
# The data directory will be specified by the command line argument if present
if len(sys.argv) > 1:
  data_dir  = sys.argv[1]
else:
  data_dir  = '/home/data'
print 'data_dir = ', data_dir

OuterLoops  = 1
Jlogplot    = False
plot_dir    = data_dir
ScalarDiags = 'diagnostics.dat'
Bg_traj     = 'LS_Oloop001_Iloop000.nc'
Anal_traj   = 'LS_Oloop%03i_Iloop000.nc'% (OuterLoops+1)
Anal_inc    = 'AnalInc.nc'
GradJo_traj = 'GradJo_Oloop001_Iloop001.nc'
Obsfile_bg  = 'Obs_001_Iloop000.dat'                     # Obs and background obs
Obsfile_an  = 'Obs_%03i_Iloop000.dat'% (OuterLoops+1)    # Obs and analysis obs
Truth_dir   = data_dir + '/Obs+Truth'
Truth_traj  = 'Truth.nc'

# Which routines do we run?
PlotBg_traj = True  # Background trajectory
PlotAn_traj = True  # Analysis trajectory
PlotAn_inc  = True   # Analysis increment
PlotJgrad   = True  # Grad of Jo trajectory
PlotTr_traj = True  # Truth trajectory
PlotBE_traj = True  # Background error trajectory
PlotAE_traj = True  # Analysis error trajectory
PlotIts     = True   # Iteration-dependent quantities
PlotBal_traj= True  # Balance characteristic trajectories (including scale-dependences)
PlotObs     = True  # Plot various observation - state histograms

# Plot all times "AllTimes" or just first and last times "FirstLast" in above options?
# The option "None" is also used to turn-off plotting altogether (and just collect data).
TimeOutput  = "FirstLast"

# Set this to output only a compact (mean error and rms error) time sequence for this da cycle
# If set to true, many of above Plotxxx and TimeOutput settings will be over-ridden.
Output_only_compact_data = True

# Set the domain dimensions
nlongs      = 360
nlevs       = 60
C_param     = 10000.0
f           = 0.0001

# Set number of bins for o-b, o-a, etc diagnostics
o_minus_model_nbins = 25

# Output type - not yet implemented - will always output to web
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'



# Over-ride settings if necessary
if Output_only_compact_data:
  PlotBg_traj = False
  PlotAn_traj = False
  PlotAn_inc  = False
  PlotJgrad   = False
  PlotTr_traj = False
  PlotBE_traj = True
  PlotAE_traj = True
  PlotIts     = True
  PlotBal_traj= False
  PlotObs     = True
  TimeOutput  = "None"


os.system('mkdir -p ' + plot_dir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
else:
  filesuffix = '.eps'


if (output_type == 'web'):
  # Set-up the html frames file
  html_file = open (plot_dir + '/Frames.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<FRAMESET ROWS=50%,50%>\n')
  html_file.write ('<FRAME src=Plots.html></FRAME>\n')
  html_file.write ('<FRAME src=Plots.html></FRAME>\n')
  html_file.write ('</FRAMESET>\n')
  html_file.write ('</html>')
  html_file.close

  # Set-up the html file
  html_file = open (plot_dir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Assimilation output</h1>\n')
  html_file.write (data_dir)
else:
  html_file = 0


# 1. Read and plot the background trajectory
# ------------------------------------------
if PlotBg_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Background trajectory</h2>\n')
  plotfields (data_dir + '/' + Bg_traj, 'Bg', output_type, html_file, plot_dir, TimeOutput)

# 2. Read and plot the analysis trajectory
# ----------------------------------------
if PlotAn_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Analysis trajectory</h2>\n')
  plotfields (data_dir + '/' + Anal_traj, 'Anal', output_type, html_file, plot_dir, TimeOutput)

# 3. Read and plot the analysis increment
# ---------------------------------------
if PlotAn_inc:
  if (output_type == 'web'):
    html_file.write ('<h2>Analysis increment</h2>\n')
  plotfields (data_dir + '/' + Anal_inc, 'AnalInc', output_type, html_file, plot_dir, TimeOutput)

# 4. Read and plot the propagation of the gradient of Jo backwards
# ----------------------------------------------------------------
if PlotJgrad:
  if (output_type == 'web'):
    html_file.write ('<h2>Gradient of Jo trajectory</h2>\n')
  plotfields (data_dir + '/' + GradJo_traj, 'GradJo', output_type, html_file, plot_dir, TimeOutput)

# 5. Read and plot the truth trajectory
# ------------------------------------------
if PlotTr_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Truth trajectory</h2>\n')
  plotfields (Truth_dir + '/' + Truth_traj, 'Tr', output_type, html_file, plot_dir, TimeOutput)

# 6. Plot the background error trajectory fields
# ----------------------------------------------
# Also return the mean background error, and RMS background error for all quantities
if PlotBE_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Background error trajectory</h2>\n')
  times, QuantityNames, TRUTH_mean, BG_mean, BG_err, TRUTH_RMS, BG_RMS, BG_RMSE = \
        plot_diff_fields (Truth_dir + '/' + Truth_traj,
                          data_dir + '/' + Bg_traj,
                          'Bg_err', output_type, html_file, plot_dir,
                          C_param, f, TimeOutput)

  # Dump the time sequence of error data to a file inside the data directory
  print 'Dumping data to file'
  dump_scalar_time_seq (times, QuantityNames, 'bg', TRUTH_mean, BG_mean, BG_err, TRUTH_RMS, BG_RMS, BG_RMSE, data_dir)
  print 'Done'

  # Plot the scalar means
  if TimeOutput  != "None":
    print 'Plotting ...'
    plot_scalar_time_seq (times, [0, times[-1]], QuantityNames, 'bg', TRUTH_mean, BG_mean, BG_err, TRUTH_RMS, BG_RMS, BG_RMSE, output_type, html_file, plot_dir)
    print 'Done'


# 7. Plot the analysis error trajectory fields
# --------------------------------------------
# Also return the mean analysis error, and RMS analysis error for all quantities
if PlotAE_traj:
  if (output_type == 'web'):
    html_file.write ('<h2>Analysis error trajectory</h2>\n')
  times, QuantityNames, TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE = \
        plot_diff_fields (Truth_dir + '/' + Truth_traj,
                          data_dir + '/' + Anal_traj,
                          'Anal_err', output_type, html_file, plot_dir,
                          C_param, f, TimeOutput)

  # Dump the time sequence of error data to a file inside the data directory
  dump_scalar_time_seq (times, QuantityNames, 'anal', TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE, data_dir)

  # Plot the scalar means
  if TimeOutput  != "None":
    plot_scalar_time_seq (times, [0, times[-1]], QuantityNames, 'anal', TRUTH_mean, ANAL_mean, ANAL_err, TRUTH_RMS, ANAL_RMS, ANAL_RMSE, output_type, html_file, plot_dir)




# Deal with the scalar diagnostics that change with iteration
# ===========================================================
if PlotIts:
  filediags = open (data_dir + '/' + ScalarDiags, 'r')
  line = filediags.readline()
  line = filediags.readline()
  innerloop = []
  outerloop = []
  totalit   = []
  Jb        = []
  Jo        = []
  J         = []
  grad      = []
  KE        = []
  BE        = []
  EE        = []
  TE        = []

  tit       = 0
  line      = filediags.readline()
  while (line != ''):
    temp = line.split()
    outerloop.append (int(temp[0]))
    innerloop.append (int(temp[1]))
    tit += 1
    totalit.append(tit)
    Jb.append (float(temp[2]))
    Jo.append (float(temp[3]))
    J.append (float(temp[4]))
    grad.append (float(temp[5]))
    if (int(temp[1]) == 0):
      # Energy output only at the start of the outer loop
      KE.append (float(temp[6]))
      BE.append (float(temp[7]))
      EE.append (float(temp[8]))
      TE.append (float(temp[9]))
    line = filediags.readline()

  filediags.close

  if (Jlogplot):
    # Replace all zeros with nans (so can do log plot)
    for element in range(len(J)):
      if (J[element] == 0.0):
        J[element] = np.nan
      if (Jb[element] == 0.0):
        Jb[element] = np.nan
      if (Jo[element] == 0.0):
        Jo[element] = np.nan

  # Find loop boundaries (to mark position of outer loop changes)
  loopboundsx, loopboundsy = OuterLoopBounds (J, Jb, Jo, outerloop, totalit)


  # 8. Plot cost function with iteration
  # ------------------------------------
  fig, ax = plt.subplots()
  ax.set_xlabel('total iteration')
  ax.set_ylabel('cost function')

  # Remove the last item as it's usually nonsense
  totalit = totalit[:-1]
  Jb      = Jb[:-1]
  Jo      = Jo[:-1]
  J       = J[:-1]
  grad    = grad[:-1]

  if (Jlogplot):
    ax.set_yscale('log')
  ax.plot(totalit[:], Jb[:], linewidth=2, color='blue', label='Jb')
  ax.plot(totalit[:], Jo[:], linewidth=2, color='red', label='Jo')
  ax.plot(totalit[:], J[:], linewidth=2, color='black', label='J')

  for loopbound in range(len(loopboundsx)):
    ax.plot(loopboundsx[loopbound][:], loopboundsy[loopbound][:], linewidth=1, color='yellow')

  #plt.scatter(iterations[:], cc[:], color='green')
  plt.title('Cost function components with iteration')
  ax.legend(loc='upper right')

  graphics_file_name = 'Plots/J' + filesuffix

  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # 9. Plot energy with iteration
  # -----------------------------
  fig, ax = plt.subplots()

  ax.set_xlabel('total iteration')
  ax.set_ylabel('energy')

  outerloopiteration = range(len(TE))
  ax.plot(outerloopiteration[:], KE[:], linewidth=2, color='red', label='Kinetic')
  ax.plot(outerloopiteration[:], BE[:], linewidth=2, color='green', label='Buoyant')
  ax.plot(outerloopiteration[:], EE[:], linewidth=2, color='blue', label='Elastic')
  ax.plot(outerloopiteration[:], TE[:], linewidth=2, color='black', label='Total')


  #plt.scatter(iterations[:], cc[:], color='green')
  plt.title('Energy at the start of each outer loop')
  ax.legend(loc='upper right')

  graphics_file_name = 'Plots/Energy' + filesuffix

  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # 10. Plot gradient of cost function with iteration
  # ------------------------------------------------
  fig, ax = plt.subplots()
  ax.set_xlabel('total iteration')
  ax.set_ylabel('gradient norm')

  ax.set_yscale('log')
  ax.plot(totalit[:], grad[:], linewidth=2, color='blue')

  #for loopbound in range(len(loopboundsx)):
  #  ax.plot(loopboundsx[loopbound][:], loopboundsy[loopbound][:], linewidth=1, color='yellow')

  #plt.scatter(iterations[:], cc[:], color='green')
  plt.title('Gradient norm with iteration')

  graphics_file_name = 'Plots/Grad' + filesuffix

  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



'''
# 11. Plot the balance characteristics with time
# ---------------------------------------------------
if PlotBal_traj:
  # Set directories
  bal_files = []                                 ; label = []
  bal_files.append(data_dir + '/' + Bg_traj)     ; label.append('Background')
  bal_files.append(data_dir + '/' + Anal_traj)   ; label.append('Analysis')
  bal_files.append(Truth_dir + '/' + Truth_traj) ; label.append('Truth')

  # Find out how many times are contained in the files
  nc_file   = Dataset(bal_files[0])
  times     = nc_file.variables['time'][:] / 3600.0
  tstep_max = nc_file.dimensions['time'].size
  nc_file.close

  Lx        = 1.5 * float(nlongs)

  # Wavenumbers and scales
  hori_wns      = np.linspace(0,nlongs-1,nlongs)
  hori_lens     = np.zeros(nlongs)
  hori_lens[1:nlongs/2+1] = Lx / hori_wns[1:nlongs/2+1]
  hori_lens[0]  = Lx
  for l in range(1, nlongs/2):
    hori_lens[nlongs-l] = hori_lens[l]

  # Define the scales of interest (km)
  scalelim = []                  ; lineone = []
  scalelim.append(100.0)         ; lineone.append('solid')
  scalelim.append(10.0)          ; lineone.append('dashed')
  scalelim.append(1.0)           ; lineone.append('dotted')


  # Set-up the arrays to store the results
  print 'Number of directories ', len(bal_files)
  print 'Number of scales      ', len(scalelim)
  gimbalresults = np.zeros((len(bal_files), tstep_max, len(scalelim)))
  himbalresults = np.zeros((len(bal_files), tstep_max, len(scalelim)))

  for fileno in range(len(bal_files)):

    bal_file = bal_files[fileno]
    print 'INPUT FILE  : ', bal_file

    for time in range(tstep_max):
      print '  Dealing with timestep ', time

      gimbal_rms, himbal_rms = balance_scale (bal_file, '', time, hori_lens, scalelim, C_param, f, nlongs, nlevs)

      # Store the diagnostics
      gimbalresults[fileno,time,:] = gimbal_rms[:]
      himbalresults[fileno,time,:] = himbal_rms[:]

  # Plot the results - separate plots for each experiment
  for fileno in range(len(bal_files)):

    html_file.write ('<h4>' + label[fileno] + '<h4>\n')

    # Plot the geostrophic imbalance
    plotfile = 'Plots/' + label[fileno] + '_Imbal' + filesuffix
    print 'Plotting to file :', plotfile

    fig, ax1 = plt.subplots()
    ax2      = ax1.twinx()

    print 'Plotting geostrophic imbalance'
    #ax1.set_ylim([0.0,0.9])
    for scale in range (0,len(scalelim)):
      ax1.plot(times, gimbalresults[fileno,:,scale], color='blue', linewidth='2', ls=lineone[scale], label='>' + str(scalelim[scale]) + 'km')
    ax1.legend(loc='upper left')
    #ax1.legend(loc=(0.3, legypos[fileno]))
    ax1.set_xlabel('time (hours)', color='black', fontsize=16)
    ax1.set_ylabel('Geostrophic imbalance', color='blue', fontsize=16)
    for tl in ax1.get_yticklabels():
      tl.set_color('blue')

    print 'Plotting hydrostatic imbalance'
    #ax2.set_ylim([0.0,0.12])
    for scale in range (0,len(scalelim)):
      ax2.plot(times, himbalresults[fileno,:,scale], color='red', linewidth='2', ls=lineone[scale], label='>' + str(scalelim[scale]) + 'km')
    ax2.legend(loc='lower right')
    #ax2.legend(loc=(0.6, legypos[fileno]))
    ax2.set_ylabel('Hydrostatic imbalance', color='red', fontsize=16)
    for tl in ax2.get_yticklabels():
      tl.set_color('red')

    plt.title(label[fileno] + ' geo and hydro imbalance', color='black', fontsize=16)
    #plt.show()
    plt.savefig(plot_dir + '/' + plotfile, bbox_inches='tight')
    plt.close()

    if (output_type == 'web'):
      html_file.write ('<br><img src=' + plotfile + ' width=300>\n')
'''



# 12. Observations, background observations, analysis observations
# ================================================================
if PlotObs:
  # Read-in the observations, background observations, and analysis observations (o = observations, bo = background model observations, to = true observations, ao = analysis model observations)
  u_o, v_o, w_o, r_o, b_o, tr_o, horwind_o, totwind_o, u_bo, v_bo, w_bo, r_bo, b_bo, tr_bo, horwind_bo, totwind_bo, u_to, v_to, w_to, r_to, b_to, tr_to, horwind_to, totwind_to = ReadObs_file (data_dir + '/' + Obsfile_bg)
  u_o, v_o, w_o, r_o, b_o, tr_o, horwind_o, totwind_o, u_ao, v_ao, w_ao, r_ao, b_ao, tr_ao, horwind_ao, totwind_ao, u_to, v_to, w_to, r_to, b_to, tr_to, horwind_to, totwind_to = ReadObs_file (data_dir + '/' + Obsfile_an)

  # Generate histograms
  if len(u_o) > 0:
    make_histogram_ob (o_minus_model_nbins, u_o, u_bo, u_ao, u_to, 'u', output_type, plot_dir, html_file)
  if len(v_o) > 0:
    make_histogram_ob (o_minus_model_nbins, v_o, v_bo, v_ao, v_to, 'v', output_type, plot_dir, html_file)
  if len(w_o) > 0:
    make_histogram_ob (o_minus_model_nbins, w_o, w_bo, w_ao, w_to, 'w', output_type, plot_dir, html_file)
  if len(r_o) > 0:
    make_histogram_ob (o_minus_model_nbins, r_o, r_bo, r_ao, r_to, 'r', output_type, plot_dir, html_file)
  if len(b_o) > 0:
    make_histogram_ob (o_minus_model_nbins, b_o, b_bo, b_ao, b_to, 'b', output_type, plot_dir, html_file)
  if len(tr_o) > 0:
    make_histogram_ob (o_minus_model_nbins, tr_o, tr_bo, tr_ao, tr_to, 'tracer', output_type, plot_dir, html_file)
  if len(horwind_o) > 0:
    make_histogram_ob (o_minus_model_nbins, horwind_o, horwind_bo, horwind_ao, horwind_to, 'horizwind', output_type, plot_dir, html_file)
  if len(totwind_o) > 0:
    make_histogram_ob (o_minus_model_nbins, totwind_o, totwind_bo, totwind_ao, totwind_to, 'totalwind', output_type, plot_dir, html_file)



if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close
