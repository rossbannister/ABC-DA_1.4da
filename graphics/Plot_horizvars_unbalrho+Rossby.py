#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to
#  1. read horizontal variances for unbalanced pressure from CVT file
#  2. read vertical lengthscales and Rossby radii from data file
#  3. sort CVT data into height order
#  4. plot horizontal variances as a function of horizontal and vertical lengthscales
#  5. plot Rossy radii on the same plot
#
# For each control variable:
#   Standard deviations
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# Ross Bannister, Dec 2016, Mar 2018
# -------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import os

# Set the input CVT file
CVT_file     = '/home/ross/DataAssim/RuthsModel/ABC_vn1.4da/Investigations/Orig_transform_order_smoothed_sigma/Exp+GB+HB-AB/MasterCalibration/CVT.nc'

# Set the intput data file (vertical lengthscales, Rossby radii of vertical modes)
lengths_file = '/home/ross/Papers/ABC_DA/HorzSpectra_processed/VerticalModeCharacteristics.dat'

# Set the domain dimensions
nlongs   = 360
nlevs    = 60

# Output directory for plots
plotdir = '.'



# Open the CVT file for reading
nc_file = Dataset(CVT_file)

# Deal with axes
cvt_order        = nc_file.variables['cvt_order'][0]
cvt_vert_opt_sym = nc_file.variables['cvt_vert_opt_sym'][0]
waven            = nc_file.variables['wavenumber'][:]
vertmode         = nc_file.variables['vert_mode'][:]
level            = nc_file.variables['level'][:] / 1000.0


# Find the max and min of the axis values
minwn = np.min(np.abs(waven))
maxwn = np.max(np.abs(waven))
if ((cvt_order == 1) and (cvt_vert_opt_sym == 1)):
  miny = np.min(np.abs(vertmode))
  maxy = np.max(np.abs(vertmode))
else:
  miny = np.min(np.abs(level))
  maxy = np.max(np.abs(level))
extent = [minwn, maxwn, miny, maxy]


# --------------------------------------------------------------
# --- Read the horizontal eigenvalues (parameter 3)          ---
# --------------------------------------------------------------
print 'Read the horizontal eigenvalues parameter 3'
evals = nc_file.variables['horizEV3'][:,:]
nc_file.close

# --------------------------------------------------------------
# --- Read the data about the vertical modes (vert length and Rossby radius) ---
# --------------------------------------------------------------
print 'Read the data about the vertical modes'

vert_length = []
Rossby_rad  = []
input_file  = open (lengths_file, 'r')
line        = input_file.readline()
for l in range(nlevs):
  line  = input_file.readline()
  split = line.split()
  vert_length.append(float(split[1]) / 1000.0)
  Rossby_rad.append(float(split[3]) / 1000.0)
input_file.close()

#for m in range(nlevs):
#  print m, vert_length[m], Rossby_rad[m]
#print
#print
#print


# --------------------------------------------------------------
# --- Sort the vertical modes into vertical length order     ---
# --------------------------------------------------------------
print 'Sorting the vertical modes'

# Set initial arrangement (order will refer to the mode indices in ascending order of vertical length)
order = range(nlevs)

# Do selection sort
for m in range(nlevs-1):
  # Find smallest of remaining values
  low_value = 9999999999.0
  low_index = m
  for mm in range(m,nlevs):
    ref = order[mm]
    if (vert_length[ref] < low_value):
      low_value = vert_length[ref]
      low_index = mm
  #print '-------------------------------------------'
  #print m, ' Lowest value = ', low_value
  if (low_index != m):
    # Need to do a swap
    ref              = order[low_index]
    order[low_index] = order[m]
    order[m]         = ref
  #print order

# Check that the lengths are ordered correctly
#for m in range(nlevs):
#  print m, vert_length[m], order[m], vert_length[order[m]]


# Rearrange data according to correct height order
evals_ordered    = np.zeros((nlevs,nlongs/2+1))
vert_len_ordered = np.zeros(nlevs)
Rossby_ordered   = np.zeros(nlevs)

for m in range(nlevs):
  ref = order[m]
  for k in range(nlongs/2+1):
    evals_ordered[m][k] = evals[ref][k]
  vert_len_ordered[m] = vert_length[ref]
  Rossby_ordered[m]   = Rossby_rad[ref]

# Check that the lengths are ordered correctly
#for m in range(nlevs):
#  print m, vert_length[m], order[m], vert_length[order[m]], vert_len_ordered[m]


# Determing the lengthscales from the horizontal wavenumbers
L         = 1.5 * float(nlongs)
horiz_len = []
horiz_len.append(1.*L)
for k in range(1,nlongs/2+1):
  horiz_len.append(L/float(k))



# --------------------------------------------------------------
# --- Plot the data                                          ---
# --------------------------------------------------------------
print 'Plotting vertical lengthscales as a function of vertical mode'
fig, ax = plt.subplots()
fig.set_size_inches(2.5, 2.0)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8) 
ax.set_xlabel('Vertical lengthscale (km)', fontsize=8)
ax.set_ylabel('Vertical mode index', fontsize=8)
plt.title('Vertical lengthscales of rho\' modes', fontsize=8)
ax.plot(vert_length[:], vertmode[:], linewidth=1, color='black')
plt.savefig('vertlengths.eps', bbox_inches='tight')
plt.close('all')


print 'Plotting horizontal spectra and Rossby radii'
# Plot filled contours of a field
fig, ax = plt.subplots()
fig.set_size_inches(5.0, 2.0)
matplotlib.rc('xtick', labelsize=8)
matplotlib.rc('ytick', labelsize=8)
cmap = cm.get_cmap('YlGnBu', 8)
ax.set_xscale('log')
cax  = ax.contourf(horiz_len[:], vert_len_ordered[:], evals_ordered[:,:], cmap=cmap)
cbar = fig.colorbar(cax, orientation='vertical', cmap=cmap)
#cax  = ax.contour(horiz_len[:], vert_len_ordered[:], evals_ordered[:,:], colors='k')
ax.plot(Rossby_ordered[:], vert_len_ordered[:], linewidth=1, color='black')
# Labels etc
ax.set_title('Sqrt horiz evals for (unbalanced) rho\'', fontsize=8)
ax.set_xlabel('horizontal lengthscale (km)', fontsize=8)
ax.set_ylabel('vertical lengthscale (km)', fontsize=8)

plt.savefig('HorizSqrtEvals_unbalrho_ordered.eps', bbox_inches='tight')
plt.close('all')

