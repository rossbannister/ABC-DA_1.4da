#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities
#  Vertical covariance matrix for balanced r with itself
#  Vertical covariance matrix for total r with balanced r
#  Vertical regression matrix used to compute regressed balanced r
#
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# Ross Bannister, Mar 2018
# -------------------------------------------------------------------


# ===================================================================

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os

# Set directory of the CVT file
datadirCVT = '/home/data'
CVT_file   = 'CVT.nc'

# Set the domain dimensions
nlevs    = 60

# Output directory for plots
plotdir = datadirCVT + '/Master_Calibration_stage3'

# Output type
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'eps'

os.system('mkdir -p ' + plotdir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
  html_file = open (plotdir + '/Plots_Regression.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Vertical regression matrices for ABC model</h1>\n')
else:
  filesuffix = '.eps'


# Open the CVT file for reading
nc_file = Dataset(datadirCVT + '/' + CVT_file)

cov_rbalrbal = nc_file.variables['cov_rbalrbal'][:,:]
cov_rtotrbal = nc_file.variables['cov_rtotrbal'][:,:]
regress_gb   = nc_file.variables['regress_gb'][:,:]
levels       = nc_file.variables['level'][:] / 1000.0

nc_file.close


# Find the max and min of the axis values
minlev = 1
maxlev = nlevs

extent = [minlev, maxlev, maxlev, minlev]

print extent

# ======================================================================
# --- Plot the vertical covariance matrix for balanced r with itself ---

if (output_type == 'web'):
  html_file.write ('<h2>Vertical error covariance matrix for balanced r with itself\n')

# Find the max and min of the matrix values
minmatrix     = np.abs(np.min(cov_rbalrbal))
maxmatrix     = np.abs(np.max(cov_rbalrbal))
maxmatrix_use = max(minmatrix,maxmatrix)
minmatrix_use = -maxmatrix_use

fig, ax = plt.subplots(1, 1, figsize=(12, 6))
cax     = ax.imshow(cov_rbalrbal, interpolation='None', cmap='seismic',
                    extent=extent, vmin=minmatrix_use, vmax=maxmatrix_use)
cbar    = fig.colorbar(cax, orientation='vertical')
ax.set_title('rbal-rbal covariance')
ax.set_xlabel('rbal level', fontsize=12)
ax.set_ylabel('rbal level', fontsize=12)

#plt.show()
plt.savefig(plotdir + '/Plots/cov_rbalrbal' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/cov_rbalrbal' + filesuffix + ' width=500></td>\n')


# ======================================================================
# --- Plot the vertical covariance matrix for total r with balanced r ---

if (output_type == 'web'):
  html_file.write ('<h2>Vertical error covariance matrix for total r with balanced r\n')

# Find the max and min of the matrix values
minmatrix     = np.abs(np.min(cov_rtotrbal))
maxmatrix     = np.abs(np.max(cov_rtotrbal))
maxmatrix_use = max(minmatrix,maxmatrix)
minmatrix_use = -maxmatrix_use

fig, ax = plt.subplots(1, 1, figsize=(12, 6))
cax     = ax.imshow(cov_rtotrbal, interpolation='None', cmap='seismic',
                    extent=extent, vmin=minmatrix_use, vmax=maxmatrix_use)
cbar    = fig.colorbar(cax, orientation='vertical')
ax.set_title('rtot-rbal covariance')
ax.set_xlabel('rbal level', fontsize=12)
ax.set_ylabel('rtot level', fontsize=12)

#plt.show()
plt.savefig(plotdir + '/Plots/cov_rtotrbal' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/cov_rtotrbal' + filesuffix + ' width=500></td>\n')


# ======================================================================
# --- Plot the vertical regression matrix used to compute regressed balanced r ---

if (output_type == 'web'):
  html_file.write ('<h2>Vertical regression matrix used to compute regressed balanced r\n')

# Find the max and min of the matrix values
minmatrix     = np.abs(np.min(regress_gb))
maxmatrix     = np.abs(np.max(regress_gb))
maxmatrix_use = max(minmatrix,maxmatrix)
minmatrix_use = -maxmatrix_use

fig, ax = plt.subplots(1, 1, figsize=(12, 6))
cax     = ax.imshow(regress_gb, interpolation='None', cmap='seismic',
                    extent=extent, vmin=minmatrix_use, vmax=maxmatrix_use)
cbar    = fig.colorbar(cax, orientation='vertical')
ax.set_title('Vertical geostrophic regression matrix')
ax.set_xlabel('level', fontsize=12)
ax.set_ylabel('level', fontsize=12)

#plt.show()
plt.savefig(plotdir + '/Plots/regress_gb' + filesuffix, bbox_inches='tight')
plt.close('all')

if (output_type == 'web'):
  html_file.write ('<td><img src=Plots/regress_gb' + filesuffix + ' width=500></td>\n')

if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close()
  print 'An html file has been created to view the figures'
  print 'Please view the following html file with your browser'
  print plotdir + '/Plots/Plots_Regression.html'
  print
