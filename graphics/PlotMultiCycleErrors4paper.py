#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read error diagnostics from each cycle, splice them, and then plot
#
# Ross Bannister, August 2018
# -------------------------------------------------------------------


import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
import os
import sys

#Base_dir  = '/home/ross/DataAssim/RuthsModel/ABC_vn1.4da/Investigations/Orig_transform_order_smoothed_sigma_noForceCalib/Exp+GB+HB-AB+Vreg/DA_cycling_Obs_rp'
Base_dir  = '/home/ross/DataAssim/RuthsModel/ABC_vn1.4da/Investigations/Orig_transform_order_smoothed_sigma_noForceCalib/Exp+GB+HB-AB-Vreg/DA_cycling_Obs_rp'
#Base_dir  = '/home/ross/DataAssim/RuthsModel/ABC_vn1.4da/Investigations/Orig_transform_order_smoothed_sigma_noForceCalib/Exp-GB+HB-AB/DA_cycling_Obs_rp'

# The file containing the list of sequential DA runs performed by the script
Exp_list    = 'ExpList.dat'

# The location of the plots
plot_dir    = Base_dir + '/ExtraPlots4Paper'

# Also plot data from the free backgound run?
Plot_freebg = True

# The full path of the file containing the forecast from the first background (relevant only if Plot_freebg = True)
FullBg_file = Base_dir + '/Master_RunNLModel_Fullbg/BgFc.nc'

# Number of output times per cycle
Noutputtimes = 6


# If requested, open the free backgound run netcdf file
if Plot_freebg:
  nc_file_bg = Dataset(FullBg_file)

# Read-in the list of files
input_file = open (Base_dir + '/' + Exp_list, 'r')
cycle_list = []
line = input_file.readline()
line = input_file.readline()
while (line != ''):
  cycle_list.append(line[:-1])  # Remove newline character
  line = input_file.readline()
input_file.close()
print 'There are ', len(cycle_list), ' cycles in this DA run'


print '  Reading from files'

# The master time sequencies for all cycles concatenated together
truth_mean      = []
data_mean       = []
freebg_mean     = []
data_mean_err   = []
freebg_mean_err = []
truth_rms       = []
data_rms        = []
data_rms_err    = []
cycle_count     = -1

for cycle_dir in cycle_list:
  cycle_count += 1
  print '  Dealing with directory : ', cycle_dir
  input_file = open (cycle_dir + '/anal.dat', 'r')
  line = input_file.readline()
  line = input_file.readline()
  # Get the times
  line = input_file.readline()
  line = input_file.readline()
  times_cycle = map(float, line.split())
  ntimes_eachcycle = len(times_cycle)

  # Get quantities that have mean and rms values and errors
  # Put into a generic structure first to discover what is available
  quantities          = []
  truth_mean_cycle    = []
  data_mean_cycle     = []
  data_mean_err_cycle = []
  truth_rms_cycle     = []
  data_rms_cycle      = []
  data_rms_err_cycle  = []
  for quantity in range(12):
    line = input_file.readline()
    quantities.append(line[:-1])  # -1 to remove newline
    line = input_file.readline()
    truth_mean_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_mean_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_mean_err_cycle.append(map(float, line.split()))
    line = input_file.readline()
    truth_rms_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_rms_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_rms_err_cycle.append(map(float, line.split()))

  # Get balance quantities that are to do with scale bands
  # These have only rms values (and we don't yet know how many bands there are)
  line = input_file.readline()
  while (line != ''):
    quantities.append(line[:-1])  # -1 to remove newline
    line = input_file.readline()
    truth_rms_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_rms_cycle.append(map(float, line.split()))
    line = input_file.readline()
    data_rms_err_cycle.append(map(float, line.split()))
    line = input_file.readline()
  Nscales = (len(data_rms_cycle) - len(data_mean_cycle)) / 2
  print '  There are ', Nscales, ' scale bands in the file'

  input_file.close()

  # Get data for the free background run if requested
  if Plot_freebg:
    freebg_mean_cycle     = []
    freebg_mean_err_cycle = []
    # Truth file
    truth_file = cycle_dir + '/Obs+Truth/Truth.nc'
    nc_file_tr = Dataset(truth_file)

    # Read and compute for each quantity
    for quantity in range(6):
      freebg  = nc_file_bg.variables[quantities[quantity]][cycle_count,:,:]
      truth   = nc_file_tr.variables[quantities[quantity]][0,:,:]
      # RMS of free background at this time
      freebg_mean_cycle.append(np.sqrt(np.mean(freebg * freebg)))
      # RMS of free background error at this time
      error = freebg - truth
      freebg_mean_err_cycle.append(np.sqrt(np.mean(error*error)))

    nc_file_tr.close


  # Create the master structures if this is the first cycle through
  # This creates lists of empty time sequences
  if len(data_mean) == 0:
    for quantity in range(len(data_mean_cycle)):
      truth_mean.append([])
      data_mean.append([])
      data_mean_err.append([])
  if len(data_rms) == 0:
    for quantity in range(len(data_rms_cycle)):
      truth_rms.append([])
      data_rms.append([])
      data_rms_err.append([])
    if Plot_freebg:
      for quantity in range(6):
        freebg_mean.append([])
        freebg_mean_err.append([])

  # Append these data to the master time sequencies
  for quantity in range(len(data_mean_cycle)):
    for time in truth_mean_cycle[quantity]:
      truth_mean[quantity].append(time)
    for time in data_mean_cycle[quantity]:
      data_mean[quantity].append(time)
    for time in data_mean_err_cycle[quantity]:
      data_mean_err[quantity].append(time)
  for quantity in range(len(data_rms_cycle)):
    for time in truth_rms_cycle[quantity]:
      truth_rms[quantity].append(time)
    for time in data_rms_cycle[quantity]:
      data_rms[quantity].append(time)
    for time in data_rms_err_cycle[quantity]:
      data_rms_err[quantity].append(time)
  if Plot_freebg:
    for quantity in range(6):
      freebg_mean[quantity].append(freebg_mean_cycle[quantity])
      freebg_mean_err[quantity].append(freebg_mean_err_cycle[quantity])



  print 'Storing the analysis data'
  data_mean_err_anal = np.asarray(data_mean_err)
  data_rms_err_anal  = np.asarray(data_rms_err)

  # Set-up time data - remember that the end time of one cycle is the start time of the next
  times = []
  dt    = (times_cycle[1]  - times_cycle[0]) / 3600.0
  Dt    = dt * float(ntimes_eachcycle - 1)
  print 'Data assimilation time step is ', dt
  print 'Cycle length is                ', Dt
  n_cycles = len(cycle_list)
  print 'There are ', n_cycles, ' cycles'
  print 'There are ', ntimes_eachcycle, ' data assimilation steps per cycle'

  for cycle in range(n_cycles):
    for step in range(ntimes_eachcycle):
      times.append(float(cycle)*Dt + float(step)*dt)

  # Set-up the cycle boundaries
  cycle_bound_times = []
  for cycle in range(n_cycles+1):
    cycle_bound_times.append(float(cycle)*Dt)


# Close the free backgound run netcdf file
if Plot_freebg:
  nc_file_bg.close




#plot_scalar_time_seq (times, cycle_bound_times, quantities, data_type, truth_mean, data_mean, data_mean_err, truth_rms, data_rms, data_rms_err, output_type, html_file, plot_dir)
#def plot_scalar_time_seq (time, cycle_bound_times, quantities, code, TRUTH_mean, DATA_mean, DATA_err, TRUTH_RMS, DATA_RMS, DATA_RMSE, output_type, html_file, plot_dir):

# ===================================================================
# ===================================================================

# ===== For u, v, w: plot the truth and analysis =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Wind component (m/s)')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot truth

ax.plot(times[:], truth_rms[0][:],  linewidth=1, ls='solid', color='red', label='true u RMS')
ax.plot(times[:], truth_rms[1][:],  linewidth=1, ls='solid', color='blue', label='true v RMS')
ax.plot(times[:], truth_rms[2][:],  linewidth=1, ls='solid', color='green', label='true w RMS')
# Plot analysis
ax.plot(times[:], data_rms[0][:],  linewidth=1, ls='dotted', color='red', label='assim u RMS')
ax.plot(times[:], data_rms[1][:],  linewidth=1, ls='dotted', color='blue', label='assim v RMS')
ax.plot(times[:], data_rms[2][:],  linewidth=1, ls='dotted', color='green', label='assim w RMS')
# Plot free background (if requested)
print 'Length of cycle_bound_times ', len(cycle_bound_times)
print 'Length of freebg_mean       ', len(freebg_mean[0])
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean[0][:],  linewidth=1, ls='dashed', color='red', label='Free bg u RMS')
  ax.plot(cycle_bound_times[:-1], freebg_mean[1][:],  linewidth=1, ls='dashed', color='blue', label='Free bg v RMS')
  ax.plot(cycle_bound_times[:-1], freebg_mean[2][:],  linewidth=1, ls='dashed', color='green', label='Free bg w RMS')

ax.legend(loc='lower right')
plt.title('RMS wind fields')
graphics_file_name = 'wind.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')

# ===== For u, v, w: plot the errors =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Wind component errors (m/s)')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot rms err
ax.plot(times[:], data_rms_err[0][:],  linewidth=1, ls='dotted', color='red', label='assim u RMS err')
ax.plot(times[:], data_rms_err[1][:],  linewidth=1, ls='dotted', color='blue', label='assim v RMS err')
ax.plot(times[:], data_rms_err[2][:],  linewidth=1, ls='dotted', color='green', label='assim w RMS err')
## Plot mean err
#ax.plot(times[:], data_mean_err[0][:],  linewidth=1, ls='dotted', color='red', label='assim u mean err')
#ax.plot(times[:], data_mean_err[1][:],  linewidth=1, ls='dotted', color='blue', label='assim v mean err')
#ax.plot(times[:], data_mean_err[2][:],  linewidth=1, ls='dotted', color='green', label='assim w mean err')
# Plot free background (if requested)
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean_err[0][:],  linewidth=1, ls='dashed', color='red', label='Free bg u RMS err')
  ax.plot(cycle_bound_times[:-1], freebg_mean_err[1][:],  linewidth=1, ls='dashed', color='blue', label='Free bg v RMS err')
  ax.plot(cycle_bound_times[:-1], freebg_mean_err[2][:],  linewidth=1, ls='dashed', color='green', label='Free bg w RMS err')
ax.legend(loc='best') #(loc='upper right')
plt.title('Wind field errors')
graphics_file_name = 'wind_err.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')


# ===== For rho': plot the truth and analysis =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Scaled density pert')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot truth
ax.plot(times[:], truth_rms[3][:],  linewidth=1, ls='solid', color='red', label='true rho'' RMS')
# Plot analysis
ax.plot(times[:], data_rms[3][:],  linewidth=1, ls='dotted', color='red', label='assim rho'' RMS')
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean[3][:],  linewidth=1, ls='dashed', color='red', label='Free bg rho'' RMS')
ax.legend(loc='best') #(loc='upper right')
plt.title('RMS scaled density pert fields')
graphics_file_name = 'rho.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')

# ===== For rho': plot the errors =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Scaled density pert errors')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot rms err
ax.plot(times[:], data_rms_err[3][:],  linewidth=1, ls='dotted', color='red', label='assim rho'' RMS err')
## Plot mean err
#ax.plot(times[:], data_mean_err[3][:],  linewidth=1, ls='dotted', color='red', label='assim rho'' mean err')
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean_err[3][:],  linewidth=1, ls='dashed', color='red', label='Free bg rho'' RMS err')
ax.legend(loc='best') #(loc='upper right')
plt.title('Scaled density pert errors')
graphics_file_name = 'rho_err.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')


# ===== For b': plot the truth and analysis =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Buoyancy pert')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot truth
ax.plot(times[:], truth_rms[4][:],  linewidth=1, ls='solid', color='blue', label='true b'' RMS')
# Plot analysis
ax.plot(times[:], data_rms[4][:],  linewidth=1, ls='dotted', color='blue', label='assim b'' RMS')
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean[4][:],  linewidth=1, ls='dashed', color='blue', label='Free bg b'' RMS')
ax.legend(loc='best') #(loc='upper right')
plt.title('RMS buoyancy pert fields')
graphics_file_name = 'b.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')


# ===== For b': plot the errors =====
fig, ax = plt.subplots()
ax.set_xlabel('time (h)')
ax.set_ylabel('Buoyancy pert errors')
fig.set_size_inches(14.0, 7.0)
for cycle in cycle_bound_times:
  plt.axvline(x=cycle, color='yellow')
# Plot rms err
ax.plot(times[:], data_rms_err[4][:],  linewidth=1, ls='dotted', color='blue', label='assim b'' RMS err')
## Plot mean err
#ax.plot(times[:], data_mean_err[4][:],  linewidth=1, ls='dotted', color='blue', label='assim b'' mean err')
if Plot_freebg:
  ax.plot(cycle_bound_times[:-1], freebg_mean_err[4][:],  linewidth=1, ls='dashed', color='blue', label='Free bg b'' RMS err')
ax.legend(loc='best') #(loc='upper right')
plt.title('Buoyancy pert errors')
graphics_file_name = 'b_err.eps'
plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
plt.close('all')

