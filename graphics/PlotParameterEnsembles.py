#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities for ensemble members
#  v_1 (normally streamfunction)
#  v_2 (normally velocity potential)
#  v_3 (normally unbalanced r)
#  v_4 (normally unbalanced b)
#  v_5 (normally unbalanced w)
#  v_6 (normally tracer)
#
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# Ross Bannister, Dec 2016, Mar 2018
# -------------------------------------------------------------------

# ===================================================================
def max_field (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

def min_field (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin

def remove_extremes (field, low, high):
  # Set values less than low to low
  # Set values larger than high to high
  new_field = []
  for line in field:
    new_line = []
    for val in line:
      if val < low: 
        new_line.append(low)
      else:
        if val > high:
          new_line.append(high)
        else:
          new_line.append(val)
    new_field.append(new_line)
  return new_field
  
  
# ===================================================================

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os

# Set directory of the input ensemble
data_dir_ConParams = '/home/data'

# Set the domain dimensions
nlongs   = 360
nlevs    = 60

# How many ensembles?
Nens       = 1
NEnsMems   = 1
Nlats      = 60

Nitems     = NEnsMems * Nlats

# Output type
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'eps'

# Output directory for plots
plotdir = data_dir_ConParams


os.system('mkdir -p ' + plotdir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
  html_file = open (plotdir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Calibration ensemble for ABC model - parameters</h1>\n')
else:
  filesuffix = '.eps'

# Names of control parameters
cp_names = []

for ens in range(Nens):

  for mem in range(Nitems):
    # Construct filename of this parameter member
    FileNamePart = 'PertParam_' + str(ens+1).zfill(3) + '_Item' + str(mem+1).zfill(3)

    # Open the parameter ensemble member file for reading
    nc_file = Dataset(data_dir_ConParams + '/' + FileNamePart + '.nc')


    # On first member that is read, find out the user options contained in this file
    # (this information is identical for all files)
    if ((ens == 0) and (mem == 0)):
      A                 = nc_file.variables['A'][0]
      B                 = nc_file.variables['B'][0]
      C                 = nc_file.variables['C'][0]
      f                 = nc_file.variables['f'][0]
      type_of_cv        = nc_file.variables['type_of_cv'][0]
      cvt_param_opt_gb  = nc_file.variables['cvt_param_opt_gb'][0]
      cvt_param_opt_hb  = nc_file.variables['cvt_param_opt_hb'][0]
      cvt_param_opt_ab  = nc_file.variables['cvt_param_opt_ab'][0]
      cvt_param_opt_reg = nc_file.variables['cvt_param_opt_reg'][0]

      if (output_type == 'web'):
        html_file.write ('<br>A = ' + str(A) + '\n')
        html_file.write ('<br>B = ' + str(B) + '\n')
        html_file.write ('<br>C = ' + str(C) + '\n')
        html_file.write ('<br>f = ' + str(f) + '\n')
        html_file.write ('<br>' + '\n')

      if (type_of_cv != 1):
        print 'Error - expecting parameters in model space'
        exit()

      # Determine names of control parameters
      cp_names.append ('streamfunction')
      cp_names.append ('velocity potential')

      if (cvt_param_opt_gb == 1):
        if (cvt_param_opt_reg == 1):
          cp_names.append ('unbalanced r (analytical balance, with regression)')
        elif (cvt_param_opt_reg == 2):
          cp_names.append ('unbalanced r (analytical balance, no regression)')
        else:
          print 'cvt_param_opt_reg option not implemented'
          exit()
      elif (cvt_param_opt_gb == 2):
        if (cvt_param_opt_reg == 1):
          cp_names.append ('unbalanced r (statistical balance, with regression)')
        elif (cvt_param_opt_reg == 2):
          cp_names.append ('unbalanced r (statistical balance, no regression)')
        else:
          print 'cvt_param_opt_reg option not implemented'
          exit()
      elif (cvt_param_opt_gb == 3):
        cp_names.append ('full r')
      else:
        print 'cvt_param_opt_gb option not implemented'
        exit()

      if (cvt_param_opt_hb == 1):
        cp_names.append ('unbalanced b (analytical balance)')
      elif (cvt_param_opt_hb == 2):
        cp_names.append ('unbalanced b (statistical balance)')
      elif (cvt_param_opt_hb == 3):
        cp_names.append ('full b')
      else:
        print 'cvt_param_opt_hb option not implemented'
        exit()

      if (cvt_param_opt_ab == 1):
        cp_names.append ('unbalanced w')
      elif (cvt_param_opt_ab == 2):
        cp_names.append ('full w')
      else:
        print 'cvt_param_opt_ab option not implemented'
        exit()

      cp_names.append ('tracer')


      if (output_type == 'web'):
        html_file.write ('<h2>Ensemble number ' + str(ens+1).zfill(3) + '</h2>\n')
        html_file.write ('<table cols=7 border=0>\n')
        html_file.write ('<tr><td></td>')
        for param in cp_names:
          html_file.write ('<td>' + param + '</td>')
        html_file.write ('</tr>\n')


    if (output_type == 'web'):
      html_file.write ('<tr>')
      html_file.write ('<td>Member ' + str(mem+1).zfill(3) + '</td>\n')


    for param in range(6):
      # -------------------------------------------------------------------
      # For this we use the autoscale
      print '-------------------'
      field = nc_file.variables['v_'+str(param+1)][:,:]
      longs = nc_file.variables['longs'][:] / 1000.0
      level = nc_file.variables['level'][:] / 1000.0
      print cp_names[param], field.shape
      minfield = min_field(field)
      maxfield = max_field(field)
      print 'Min value of this field ', minfield
      print 'Max value of this field ', maxfield
      mn       = np.mean(field)
      err      = field - mn
      rms      = np.sqrt(np.mean(err * err))
      print 'mean of this field      ', mn
      print 'RMS of this field       ', rms
      #print 'Contour levels: ', con_levs

      matplotlib.rc('xtick', labelsize=16)
      matplotlib.rc('ytick', labelsize=16)
      fig, ax = plt.subplots()
      #cmap    = cm.get_cmap('Greys', 11)
      cmap    = cm.get_cmap('seismic', 11)
      cax     = ax.contourf(longs, level, field, cmap=cmap)
      cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
      if (rms > 0.0):
        cax     = ax.contour (longs, level, field, colors='k')
      # Labels etc
      ax.set_title(cp_names[param] + '\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
      ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
      ax.set_ylabel('Vertical distance (km)', fontsize=16)
      #plt.show()
      plt.savefig(plotdir + '/Plots/' + FileNamePart + '_v' + str(param+1) + filesuffix, bbox_inches='tight')
      plt.close('all')

      if (output_type == 'web'):
        html_file.write ('<td><img src=Plots/' + FileNamePart + '_v' + str(param+1) + filesuffix + ' width=300></td>\n')


    if (output_type == 'web'):
      html_file.write ('\n</tr>\n')

    nc_file.close


  if (output_type == 'web'):
    html_file.write ('</table>\n')


if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close()
  print 'An html file has been created to view the figures'
  print 'Please view the following html file with your browser'
  print plotdir + '/Plots.html'
  print
