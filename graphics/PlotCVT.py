#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities from the file
# describing the control variable transform (CVT)
#
# For each control variable:
#   Standard deviations
#   Vertical modes
#   Eigenvalues of vertical modes
#   Eigenvalues of horizontal modes
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

# Set directory of the input ensemble
data_dir_CVT = 'CVT_dir'
CVT_file     = 'CVT.nc'

# Set the domain dimensions
nlongs   = 360
nlevs    = 60

# Output type
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'

# Output directory for plots
plotdir = data_dir_CVT


os.system('mkdir -p ' + plotdir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
  html_file = open (plotdir + '/Plots_CVT.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Control Variable Transform Diagnostics</h1>\n')
else:
  filesuffix = '.eps'

# Names of control parameters
cp_names = []

# Open the CVT file for reading
nc_file = Dataset(data_dir_CVT + '/' + CVT_file)


# Find out the user options contained in this file
A                 = nc_file.variables['A'][0]
B                 = nc_file.variables['B'][0]
C                 = nc_file.variables['C'][0]
f                 = nc_file.variables['f'][0]
cvt_order         = nc_file.variables['cvt_order'][0]
cvt_param_opt_gb  = nc_file.variables['cvt_param_opt_gb'][0]
cvt_param_opt_hb  = nc_file.variables['cvt_param_opt_hb'][0]
cvt_param_opt_ab  = nc_file.variables['cvt_param_opt_ab'][0]
cvt_param_opt_reg = nc_file.variables['cvt_param_opt_reg'][0]
cvt_vert_opt_sym  = nc_file.variables['cvt_vert_opt_sym'][0]
cvt_stddev_opt    = nc_file.variables['cvt_stddev_opt'][0]

if (output_type == 'web'):
  html_file.write ('<br>A = ' + str(A) + '\n')
  html_file.write ('<br>B = ' + str(B) + '\n')
  html_file.write ('<br>C = ' + str(C) + '\n')
  html_file.write ('<br>f = ' + str(f) + '\n')
  html_file.write ('<br>' + '\n')


# --------------------------------------------------------------
# --- Determine names of control parameters based on options ---
# --------------------------------------------------------------
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


# --------------------------------------------------------------
# --- Plot the standard deviations                           ---
# --------------------------------------------------------------
print 'Plotting the standard deviations of the parameters'

if (output_type == 'web'):
  html_file.write ('<h2>Standard deviations for each control parameter</h2>\n')

if (cvt_stddev_opt == 1):
  # ===== Standard deviation constant for each parameter
  for param in range(6):
    stddev = nc_file.variables['v_'+str(param+1)][0]
    if (output_type == 'web'):
      html_file.write ('<br>Parameter ' + str(param+1) + ' (' + cp_names[param] +
                       '): ' + str(stddev) + '\n')
    else:
      print 'Parameter ' + str(param+1) + ' (' + cp_names[param] + '): ' + str(stddev)


elif (cvt_stddev_opt == 2):
  # ===== Standard deviation function of height for each parameter
  if (output_type == 'web'):
    html_file.write ('<table cols=6 border=0>\n')
    html_file.write ('<tr>')

  for param in range(6):
    stddev = nc_file.variables['sigma'+str(param+1)][:]
    level  = nc_file.variables['level'][:] / 1000.0

    fig, ax = plt.subplots()
    ax.set_xlabel('standard deviation', fontsize=16)
    ax.set_ylabel('height (km)', fontsize=16)
    plt.title('Std dev of ' + cp_names[param], fontsize=16)
    ax.plot(stddev[:], level[:], linewidth=2, color='black')
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plotdir + '/Plots/sigma' + str(param+1) + filesuffix)  # Comment out to not save the file
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=' + 'Plots/sigma' + str(param+1) + filesuffix + ' width=300></td>\n')

  if (output_type == 'web'):
    html_file.write ('</tr>\n')
    html_file.write ('</table>\n')


elif (cvt_stddev_opt == 3):
  # ===== Standard deviation function of longitude and height for each parameter
  # ****** THIS CODE HAS NOT BEEN TESTED YET! ******
  if (output_type == 'web'):
    html_file.write ('<table cols=6 border=0>\n')
    html_file.write ('<tr>')

  for param in range(6):
    stddev    = nc_file.variables['sigma'+str(param+1)][:]
    longitude = nc_file.variables['longs'][:] / 1000.0
    level     = nc_file.variables['level'][:] / 1000.0
    
    mn       = np.mean(stddev)
    err      = stddev - mn
    rms      = np.sqrt(np.mean(err * err))

    fig, ax = plt.subplots(1, 1, figsize=(6, 6),
                           subplot_kw={'xticks': [], 'yticks': []})
    cax     = ax.imshow(stddev, interpolation='None', cmap='YlGnBu')
    cbar    = fig.colorbar(cax, orientation='vertical')
    ax.set_title('Std dev of\n' + cp_names[param], fontsize=10)
    ax.set_xlabel('longitudinal distance (km)', fontsize=10)
    ax.set_ylabel('height (km)', fontsize=10)
    #plt.show()
    plt.savefig(plotdir + '/Plots/sigma' + str(param+1) + filesuffix)  # Comment out to not save the file
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=' + 'Plots/sigma' + str(param+1) + filesuffix + ' width=300></td>\n')

  if (output_type == 'web'):
    html_file.write ('</tr>\n')
    html_file.write ('</table>\n')

else:
  print 'cvt_stddev_opt option not implemented'
  exit()




# --------------------------------------------------------------
# --- Plot the vertical modes                                ---
# --------------------------------------------------------------
print 'Plotting the vertical modes of each parameter'
if (output_type == 'web'):
  html_file.write ('<h2>Vertical modes for each control parameter</h2>\n')
  html_file.write ('<table cols=' + str(nlevs) + ' border=0>\n')

if (cvt_order == 1):
  # Classic ordering of vertical and horizontal transforms
  level = nc_file.variables['level'][:] / 1000.0
  for param in range(6):
    print 'Parameter ', param+1
    if (output_type == 'web'):
      html_file.write ('<tr>\n')
    for vm in range(nlevs):
      print '  Vertical Mode ', vm+1
      thismode = nc_file.variables['vertmode'+str(param+1)][vm,:]
      fig, ax = plt.subplots(1, 1, figsize=(3.5, 6))
      for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                   ax.get_xticklabels() + ax.get_yticklabels()):
        item.set_fontsize(10)
      ax.set_ylabel('height (km)', fontsize=10)
      plt.title('Vert mode ' + str(vm+1) + ' for\n' + cp_names[param], fontsize=10)
      ax.plot(thismode[:], level[:], linewidth=2, color='black')
      #plt.show()             # Comment out to not display plot on the screen
      plt.savefig(plotdir + '/Plots/vertmode' + str(vm+1) + '_' + str(param+1) + filesuffix)  # Comment out to not save the file
      plt.close('all')
      if (output_type == 'web'):
        html_file.write ('<td><img src=' + 'Plots/vertmode' + str(vm+1) + '_' + str(param+1) + filesuffix + ' width=300></td>\n')
    if (output_type == 'web'):
      html_file.write ('</tr>\n')
      
elif (cvt_order == 2):
  # Reversed ordering of vertical and horizontal transforms
  # ****** THIS CODE HAS NOT BEEN TESTED YET! ******
  level = nc_file.variables['level'][:] / 1000.0
  waven = nc_file.variables['wavenumber'][:]
  for param in range(6):
    print 'Parameter ', param+1
    if (output_type == 'web'):
      html_file.write ('<tr>\n')
    for vm in range(nlevs):
      print '  Vertical Mode ', vm+1
      thismode = nc_file.variables['vertmode'+str(param+1)][vm,:,:]
      fig, ax = plt.subplots(1, 1, figsize=(6, 6),
                             subplot_kw={'xticks': [], 'yticks': []})
      cax     = ax.imshow(thismode, interpolation='None', cmap='YlGnBu')
      cbar    = fig.colorbar(cax, orientation='vertical')
      ax.set_title('Vert mode ' + str(vm+1) + ' for\n' + cp_names[param], fontsize=10)
      ax.set_xlabel('wavenumber', fontsize=10)
      ax.set_ylabel('height (km)', fontsize=10)
      #plt.show()             # Comment out to not display plot on the screen
      plt.savefig(plotdir + '/Plots/vertmode' + str(vm+1) + '_' + str(param+1) + filesuffix)  # Comment out to not save the file
      plt.close('all')
      if (output_type == 'web'):
        html_file.write ('<td><img src=' + 'Plots/vertmode' + str(vm+1) + '_' + str(param+1) + filesuffix + ' width=300></td>\n')
    if (output_type == 'web'):
      html_file.write ('</tr>\n')
else:
  print 'cvt_order option not implemented'
  exit()

if (output_type == 'web'):
  html_file.write ('</table>\n')



# --------------------------------------------------------------
# --- Plot the vertical eigenvalues                          ---
# --------------------------------------------------------------
print 'Plotting the vertical eigenvalues of each parameter'
if (output_type == 'web'):
  html_file.write ('<h2>Sqrt of vertical eigenvalues for each control parameter</h2>\n')
  html_file.write ('<table cols=6 border=0>\n')
  html_file.write ('<tr>\n')

if (cvt_order == 1):
  # Classic ordering of vertical and horizontal transforms
  vertmode = nc_file.variables['vert_mode'][:]
  if (output_type == 'web'):
    html_file.write ('<tr>\n')
  for param in range(6):
    print 'Parameter ', param+1
    evals = nc_file.variables['vertEV'+str(param+1)][:]
    fig, ax = plt.subplots(figsize=(4, 8))
    for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                 ax.get_xticklabels() + ax.get_yticklabels()):
      item.set_fontsize(10)
    plt.title('Sqrt vert evals  for\n' + cp_names[param], fontsize=10)
    ax.set_xlabel('sqrt eigenvalue', fontsize=10)
    ax.set_ylabel('vertical mode index', fontsize=10)
    ax.plot(evals[:], vertmode[:], linewidth=2, color='black')
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plotdir + '/Plots/vertevals_' + str(param+1) + filesuffix)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<td><img src=' + 'Plots/vertevals_' + str(param+1) + filesuffix + ' width=300></td>\n')
  if (output_type == 'web'):
    html_file.write ('</tr>\n')
      
elif (cvt_order == 2):
  # Reversed ordering of vertical and horizontal transforms
  # ****** THIS CODE HAS NOT BEEN TESTED YET! ******
  waven    = nc_file.variables['wavenumber'][:]
  vertmode = nc_file.variables['vert_mode'][:]
  if (output_type == 'web'):
    html_file.write ('<tr>\n')
  for param in range(6):
    print 'Parameter ', param+1
    evals = nc_file.variables['vertEV'+str(param+1)][:,:]
    fig, ax = plt.subplots(1, 1, figsize=(6, 6),
                           subplot_kw={'xticks': [], 'yticks': []})
    cax     = ax.imshow(evals, interpolation='None', cmap='YlGnBu')
    cbar    = fig.colorbar(cax, orientation='vertical')
    ax.set_title('Sqrt vert evals  for\n' + cp_names[param], fontsize=10)
    ax.set_xlabel('wavenumber', fontsize=10)
    ax.set_ylabel('vertical mode index', fontsize=10)
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plotdir + '/Plots/vertevals_' + str(param+1) + filesuffix)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<td><img src=' + 'Plots/vertevals_' + str(param+1) + filesuffix + ' width=300></td>\n')
  if (output_type == 'web'):
    html_file.write ('</tr>\n')
else:
  print 'cvt_order option not implemented'
  exit()

if (output_type == 'web'):
  html_file.write ('</tr>\n')
  html_file.write ('</table>\n')


# --------------------------------------------------------------
# --- Plot the horizontal eigenvalues                        ---
# --------------------------------------------------------------
print 'Plotting the horizontal eigenvalues of each parameter'

if (output_type == 'web'):
  html_file.write ('<h2>Sqrt of horizontal eigenvalues for each control parameter</h2>\n')
  html_file.write ('<table cols=6 border=0>\n')
  html_file.write ('<tr>\n')

waven    = nc_file.variables['wavenumber'][:]
vertmode = nc_file.variables['vert_mode'][:]
level    = nc_file.variables['level'][:] / 1000.0

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

if (output_type == 'web'):
  html_file.write ('<tr>\n')
for param in range(6):
  print 'Parameter ', param+1
  evals = nc_file.variables['horizEV'+str(param+1)][:,:]
  fig, ax = plt.subplots(1, 1, figsize=(6, 3))
  cax     = ax.imshow(evals, interpolation='None', cmap='YlGnBu', extent=extent)
  cbar    = fig.colorbar(cax, orientation='vertical')
  ax.set_title('Sqrt horiz evals  for\n' + cp_names[param], fontsize=10)
  ax.set_xlabel('wavenumber', fontsize=10)
  if ((cvt_order == 1) and (cvt_vert_opt_sym == 1)):
    ax.set_ylabel('vertical mode index', fontsize=10)
  else:
    ax.set_ylabel('height ', fontsize=10)
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plotdir + '/Plots/horizevals_' + str(param+1) + filesuffix)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + 'Plots/horizevals_' + str(param+1) + filesuffix + ' width=300></td>\n')
if (output_type == 'web'):
  html_file.write ('</tr>\n')


nc_file.close

if (output_type == 'web'):
  html_file.write ('</html>')
  html_file.close()
  print 'An html file has been created to view the figures'
  print 'Please view the following html file with your browser'
  print plotdir + '/Plots_CVT.html'
  print
