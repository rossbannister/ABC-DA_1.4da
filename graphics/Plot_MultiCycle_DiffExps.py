#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read error diagnostics from each cycle, splice them
# Do this for TWO different experiments.  Take the difference and plot.
# Also output differences to a file for possible later processing
#
# Computes the difference experiment 2 diags - experiment 1 diags
#
# Ross Bannister, September 2018
# -------------------------------------------------------------------


import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors, cm
import matplotlib
from Routines4PlotAssimDiags import plot_scalar_time_seq, dump_scalar_time_seq
import os
import sys

Core_loc  = '/home/data'
Exp1      = 'Exp+GB+HB-AB'  # This is the control experiment
Exp2      = 'Exp-GB+HB-AB'
Quantity  = 'rp'
Base_dirs = []
Base_dirs.append(Core_loc + '/' + Exp1 + '/DA_cycling_Obs_' + Quantity)
Base_dirs.append(Core_loc + '/' + Exp2 + '/DA_cycling_Obs_' + Quantity)
plot_dir  = Core_loc + '/' + Exp2 + '__minus__' + Exp1 + '__' + Quantity
Text      = 'POSITIVE VALUES INDICATE THAT THE CONTROL EXPERIMENT IS BETTER'
data_type = 'anal'   # anal or bg

print 'Experiment 1: ', Base_dirs[0]
print 'Experiment 2: ', Base_dirs[1]

# The file containing the list of sequential DA runs performed by the script
Exp_list    = 'ExpList.dat'

# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'

os.system('mkdir -p ' + plot_dir + '/Plots')
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
  html_file.write ('<h1>Difference diagnostics for two DA experiments</h1>\n')
  html_file.write (Text + '<br>\n')
  html_file.write ('CONTROL EXP : ' + Base_dirs[0] + '<br>\n')
  html_file.write ('NEW EXP : : : ' + Base_dirs[1] + '<br>\n')


cycle_lists = []

for experiment in Base_dirs:
  # Read-in the list of files from this experiment
  input_file = open (experiment + '/' + Exp_list, 'r')
  cycle_list = []
  line = input_file.readline()
  line = input_file.readline()
  while (line != ''):
    cycle_list.append(line[:-1])  # Remove newline character
    line = input_file.readline()
  input_file.close()
  print 'There are ', len(cycle_list), ' cycles in this experiment'
  cycle_lists.append(cycle_list)


if len(cycle_list[0]) != len(cycle_list[1]):
  print 'Error - number of cycles in experiments must equal.'
  exit()


# Set-up the structures to hold the data for each experiment
# The master time sequencies for all cycles concatenated together
truth_mean    = []
data_mean     = []
data_mean_err = []
truth_rms     = []
data_rms      = []
data_rms_err  = []


# Read the data for the DA cycles for each of the two experiments

for experiment in range(2):

  truth_mean_this_exp    = []
  data_mean_this_exp     = []
  data_mean_err_this_exp = []
  truth_rms_this_exp     = []
  data_rms_this_exp      = []
  data_rms_err_this_exp  = []

  for cycle_dir in cycle_lists[experiment]:
    print '  Dealing with directory : ', cycle_dir
    input_file = open (cycle_dir + '/' + data_type + '.dat', 'r')
    line = input_file.readline()
    line = input_file.readline()
    # Get the times
    line = input_file.readline()
    line = input_file.readline()
    times_cycle = map(float, line.split())
    ntimes_eachcycle = len(times_cycle)
    #print '    times_cycle:', times_cycle

    # Get quantities that have mean and rms values and errors
    # Put into a generic structure first to discover what is available
    quantities               = []
    truth_mean_this_cycle    = []
    data_mean_this_cycle     = []
    data_mean_err_this_cycle = []
    truth_rms_this_cycle     = []
    data_rms_this_cycle      = []
    data_rms_err_this_cycle  = []
    for quantity in range(12):
      line = input_file.readline()
      quantities.append(line[:-1])  # -1 to remove newline
      line = input_file.readline()
      truth_mean_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_mean_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_mean_err_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      truth_rms_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_rms_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_rms_err_this_cycle.append(map(float, line.split()))

    # Get balance quantities that are to do with scale bands
    # These have only rms values (and we don't yet know how many bands there are)
    line = input_file.readline()
    while (line != ''):
      quantities.append(line[:-1])  # -1 to remove newline
      line = input_file.readline()
      truth_rms_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_rms_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
      data_rms_err_this_cycle.append(map(float, line.split()))
      line = input_file.readline()
    Nscales = (len(data_rms_this_cycle) - len(data_mean_this_cycle)) / 2
    print '  There are ', Nscales, ' scale bands in the file'

    input_file.close()

    # Create the master structures if this is the first cycle through
    # This creates lists of empty time sequences
    if len(data_mean_this_exp) == 0:
      for quantity in range(len(data_mean_this_cycle)):
        truth_mean_this_exp.append([])
        data_mean_this_exp.append([])
        data_mean_err_this_exp.append([])
    if len(data_rms_this_exp) == 0:
      for quantity in range(len(data_rms_this_cycle)):
        truth_rms_this_exp.append([])
        data_rms_this_exp.append([])
        data_rms_err_this_exp.append([])

    # Append these data to the master time sequencies for this experiment
    for quantity in range(len(data_mean_this_cycle)):
      for time in truth_mean_this_cycle[quantity]:
        truth_mean_this_exp[quantity].append(time)
      for time in data_mean_this_cycle[quantity]:
        data_mean_this_exp[quantity].append(time)
      for time in data_mean_err_this_cycle[quantity]:
        data_mean_err_this_exp[quantity].append(time)
    for quantity in range(len(data_rms_this_cycle)):
      for time in truth_rms_this_cycle[quantity]:
        truth_rms_this_exp[quantity].append(time)
      for time in data_rms_this_cycle[quantity]:
        data_rms_this_exp[quantity].append(time)
      for time in data_rms_err_this_cycle[quantity]:
        data_rms_err_this_exp[quantity].append(time)

  # Append to the structure that contains both experiments (and convert to numpy arrays)
  truth_mean.append(np.asarray(truth_mean_this_exp))
  data_mean.append(np.asarray(data_mean_this_exp))
  data_mean_err.append(np.asarray(data_mean_err_this_exp))
  truth_rms.append(np.asarray(truth_rms_this_exp))
  data_rms.append(np.asarray(data_rms_this_exp))
  data_rms_err.append(np.asarray(data_rms_err_this_exp))


# Take the difference between each experiment
truth_mean.append(truth_mean[1] - truth_mean[0])
data_mean.append(data_mean[1] - data_mean[0])
data_mean_err.append(data_mean_err[1] - data_mean_err[0])
truth_rms.append(truth_rms[1] - truth_rms[0])
data_rms.append(data_rms[1] - data_rms[0])
data_rms_err.append(data_rms_err[1] - data_rms_err[0])


# Set-up some stuff for plotting
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
#print 'The times are'
#print times

# Set-up the cycle boundaries
cycle_bound_times = []
for cycle in range(n_cycles+1):
  cycle_bound_times.append(float(cycle)*Dt)


# Plot the difference data
print 'Plotting for ', data_type
plot_scalar_time_seq (times, cycle_bound_times, quantities, data_type, truth_mean[2], data_mean[2], data_mean_err[2], truth_rms[2], data_rms[2], data_rms_err[2], output_type, html_file, plot_dir)

# Dump the data to a file
print 'Dumping data for ', data_type
dump_scalar_time_seq (times, quantities, data_type, truth_mean[2], data_mean[2], data_mean_err[2], truth_rms[2], data_rms[2], data_rms_err[2], plot_dir)




if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close
