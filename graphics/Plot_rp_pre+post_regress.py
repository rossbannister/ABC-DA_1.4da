#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities
#  total density pert
#  balanced density pert before regression
#  balanced density pert after regression
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

# Set directory of the pert ensemble
data_dir_stage2 = '/home/Master_Calibration/Master_Calibration_stage2'
# Set directory of the example balanced r_prime fields
data_dir_stage4 = '/home/Master_Calibration/Master_Calibration_stage4'

# Set the domain dimensions
nlongs   = 360
nlevs    = 60
tstep    = 0


# Output type
# Option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'eps'

# Output directory for plots
plotdir = data_dir_stage4


os.system('mkdir -p ' + plotdir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
  html_file = open (plotdir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Scaled density perts before and after vertical regression</h1>\n')
else:
  filesuffix = '.eps'



# Read-in data
# --------------

# Total r_prime
nc_file     = Dataset(data_dir_stage2 + '/PertABC_Ens001_Item001.nc')
field_total = nc_file.variables['r_prime'][tstep,:,:]
longs_v     = nc_file.variables['longs_v'][:] / 1000.0
half_level  = nc_file.variables['half_level'][:] / 1000.0
nc_file.close


# Pre-regression balanced r_prime data
nc_file     = Dataset(data_dir_stage4 + '/r_b_preregress.nc')
field_pre   = nc_file.variables['r_b_preregress'][:,:]
nc_file.close

# Post-regression balanced r_prime data
nc_file     = Dataset(data_dir_stage4 + '/r_b_postregress.nc')
field_post  = nc_file.variables['r_b_postregress'][:,:]
nc_file.close

# Find the min and max values
minfield = min(min_field(field_total), min_field(field_pre), min_field(field_post))
maxfield = max(max_field(field_total), max_field(field_pre), max_field(field_post))
print 'Overall min ', minfield
print 'Overall max ', maxfield

con_levs = np.linspace(minfield, maxfield, 11)
#cmap    = cm.get_cmap('Greys', len(con_levs))
cmap    = cm.get_cmap('seismic', len(con_levs))

# These two lines for each field are tricks for making sure that the scales remain the same for each plot
field_total[0][0] = minfield ; field_total[0][1] = maxfield
field_total       = remove_extremes (field_total, minfield, maxfield)

field_pre[0][0]   = minfield ; field_pre[0][1]   = maxfield
field_pre         = remove_extremes (field_pre, minfield, maxfield)

field_post[0][0]  = minfield ; field_post[0][1]  = maxfield
field_post        = remove_extremes (field_post, minfield, maxfield)


# Plot the total data
# ----------------------------
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
fig, ax = plt.subplots()
cax     = ax.contourf(longs_v, half_level, field_total, clevs=con_levs, cmap=cmap)
cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
cax     = ax.contour (longs_v, half_level, field_total, clevs=con_levs, colors='k')
# Labels etc
ax.set_title('r_prime', fontsize=16)
ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
ax.set_ylabel('Vertical distance (km)', fontsize=16)
#plt.show()
plt.savefig(plotdir + '/Plots/r_total' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/r_total' + filesuffix + ' width=300></td>\n')


# Plot the pre-regression data
# ----------------------------
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
fig, ax = plt.subplots()
cax     = ax.contourf(longs_v, half_level, field_pre, clevs=con_levs, cmap=cmap)
cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
cax     = ax.contour (longs_v, half_level, field_pre, clevs=con_levs, colors='k')
# Labels etc
ax.set_title('r_b, pre-regression', fontsize=16)
ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
ax.set_ylabel('Vertical distance (km)', fontsize=16)
#plt.show()
plt.savefig(plotdir + '/Plots/r_b_pre' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/r_b_pre' + filesuffix + ' width=300></td>\n')


# Plot the post-regression data
# ----------------------------
matplotlib.rc('xtick', labelsize=16)
matplotlib.rc('ytick', labelsize=16)
fig, ax = plt.subplots()
cax     = ax.contourf(longs_v, half_level, field_post, clevs=con_levs, cmap=cmap)
cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
cax     = ax.contour (longs_v, half_level, field_post, clevs=con_levs, colors='k')
# Labels etc
ax.set_title('r_b, post-regression', fontsize=16)
ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
ax.set_ylabel('Vertical distance (km)', fontsize=16)
#plt.show()
plt.savefig(plotdir + '/Plots/r_b_post' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/r_b_post' + filesuffix + ' width=300></td>\n')



if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close()
  print 'An html file has been created to view the figures'
  print 'Please view the following html file with your browser'
  print plotdir + '/Plots.html'
  print
