#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot total energy to show loss of conservation
#
# The plotting code is repeated for data in specified data directories
# Results are shown on a single plot for all data
#
# Please edit the input details (e.g. location of data files) below.
#
# Ross Bannister, Dec 2016
# -------------------------------------------------------------------

# ===================================================================
import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

# Set directories (each is a different experiment)
data_dir = []
data_dir.append('../examples/Master_RunNLModel')
# Can add more directories in the same way

# Set legends for each directory of data
legend = []              ;  thisls = []              ;  thiscol = []
legend.append('Ref')     ;  thisls.append('solid')   ;  thiscol.append('black')
# Possibilities (e.g.): solid, dotted, ...

DoPlot = True

print thiscol
datafile = 'ABC_Diagnostics.dat'

plotdir  = '../examples/Master_RunNLModel'
os.system('mkdir -p ' + plotdir + '/Plots')

# The last time considered (seconds)
tstep    = 10800.0

if (DoPlot):
  # Set-up the plot
  fig, ax = plt.subplots()

counter = -1

for datadir in data_dir:

  counter        += 1
  input_file_name = datadir + '/' + datafile
  input_file      = open (input_file_name, 'r')

  print 'Dealing with ', input_file_name

  energy = []
  t      = 0.0
  time   = []

  while t < tstep:
    line     = input_file.readline()
    data     = line.split()
    e        = float(data[5])
    if (t == 0.0):
      e0 = e
    t        = float(data[1])
    energy.append(e/e0)
    time.append(t/3600.0)

  input_file.close()

  print '  Percentage of energy lost : ', 100.0 * (1.0 - e/e0)
  

  if (DoPlot):
    # Add this line to the plot
    cax = ax.plot(time, energy, label=legend[counter], ls=thisls[counter], color=thiscol[counter], linewidth=2)

if (DoPlot):
  # Labels etc
  ax.set_title('(b) Relative total energy', fontsize=16)
  ax.set_xlabel('time (h)', fontsize=16)
  ax.set_ylabel('Total energy (E/E0)', fontsize=16)
  ax.legend(loc='lower left', fontsize=16)
  ax.set_ylim([0.92, 1.01])
  #plt.show()
  plt.savefig(plotdir + '/Plots/Energy.eps')
  plt.close('all')
