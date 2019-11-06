#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities for ensemble members
#  Zonal wind
#  Medidional wind
#  Vertical wind speed
#  density pert (r_prime)
#  b_prime
#  tracer
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
data_dir_ABC = '/home/data'

# Set the domain dimensions
nlongs   = 360
nlevs    = 60
tstep    = 0

# How many ensembles?
Nens       = 1
NEnsMems   = 1
Nlats      = 60

Nitems     = NEnsMems * Nlats

# Are the fields to be read-in full fields or perturbations?
# Choose 'Full' for full fields, or 'Pert' for perturbations
FullOrPert = 'Full'

# Output type
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'

# Output directory for plots
plotdir = data_dir_ABC


os.system('mkdir -p ' + plotdir + '/Plots')
if (output_type == 'web'):
  filesuffix = '.png'
  html_file = open (plotdir + '/Plots.html', 'w')
  html_file.write ('<html>\n')
  html_file.write ('<h1>Calibration ensemble for ABC model</h1>\n')
  html_file.write ('Full field or perturbations: ' + FullOrPert + '\n')
else:
  filesuffix = '.eps'


for ens in range(Nens):

  if (output_type == 'web'):
    html_file.write ('<h2>Ensemble number ' + str(ens+1).zfill(3) + '</h2>\n')
    html_file.write ('<table cols=7 border=0>\n')
    html_file.write ('<tr><td>.</td><td>Zonal wind</td><td>Meridional wind</td><td>Vertical wind</td><td>r_prime</td><td>b_prime</td><td>tracer</td></tr>\n')

  for mem in range(Nitems):
    # Construct filename of this member/pert
    if (FullOrPert == 'Full'):
      FileNamePart = 'FC_Ens' + str(ens+1).zfill(3) + '_Item' + str(mem+1).zfill(3)
    else:
      FileNamePart = 'PertABC_Ens' + str(ens+1).zfill(3) + '_Item' + str(mem+1).zfill(3)


    # Open the full or perturbation ensemble member file for reading
    nc_file = Dataset(data_dir_ABC + '/' + FileNamePart + '.nc')

    if (output_type == 'web'):
      html_file.write ('<tr>')
      html_file.write ('<td>Member ' + str(mem+1).zfill(3) + '</td>\n')

    # Deal with the zonal wind
    # -------------------------------------------------------------------
    # For this we use the autoscale
    print '-------------------'
    field      = nc_file.variables['u'][tstep,:,:]
    longs_u    = nc_file.variables['longs_u'][:] / 1000.0
    half_level = nc_file.variables['half_level'][:] / 1000.0
    print ' u: ', field.shape
    #minlev = -0.04
    #maxlev =  0.04
    #con_levs = np.linspace(minlev, maxlev, 11)
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
    cax     = ax.contourf(longs_u, half_level, field, cmap=cmap)
    cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
    cax     = ax.contour (longs_u, half_level, field, colors='k')
    # Labels etc
    ax.set_title('Zonal wind\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_u' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_u' + filesuffix + ' width=300></td>\n')



    # Deal with the meridional wind
    # -------------------------------------------------------------------
    # For this we use the autoscale
    print '-------------------'
    field      = nc_file.variables['v'][tstep,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    half_level = nc_file.variables['half_level'][:] / 1000.0
    print ' v: ', field.shape
    #minlev = -0.04
    #maxlev =  0.04
    #con_levs = np.linspace(minlev, maxlev, 11)
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
    cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
    cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
    cax     = ax.contour (longs_v, half_level, field, colors='k')
    # Labels etc
    ax.set_title('meridional wind\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_v' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_v' + filesuffix + ' width=300></td>\n')



    # Deal with the vertical wind
    # -------------------------------------------------------------------
    # For this, we plot -0.4 to 0.4, but mask out anything above 1, below -1
    print '-------------------'
    field      = nc_file.variables['w'][tstep,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    full_level = nc_file.variables['full_level'][:] / 1000.0
    print 'w : ', field.shape
    minlev = -0.4
    maxlev =  0.4
    con_levs = np.linspace(minlev, maxlev, 11)
    minfield = min_field(field)
    maxfield = max_field(field)
    print 'Min value of this field ', minfield
    print 'Max value of this field ', maxfield
    mn       = np.mean(field)
    err      = field - mn
    rms      = np.sqrt(np.mean(err * err))
    print 'RMS of this field       ', rms
    print 'Contour levels: ', con_levs
    # These two lines are tricks for making sure that the scales remain the same for all plots
    field[0][0] = minlev ; field[0][1] = maxlev
    field = remove_extremes (field, minlev, maxlev)

    matplotlib.rc('xtick', labelsize=16)
    matplotlib.rc('ytick', labelsize=16)
    fig, ax = plt.subplots()
    cmap    = cm.get_cmap('seismic', len(con_levs))
    cax     = ax.contourf(longs_v, full_level, field, clevs=con_levs, cmap=cmap, extend='both')
    cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
    cax     = ax.contour (longs_v, full_level, field, clevs=con_levs, colors='k')
    # Labels etc
    ax.set_title('vertical wind\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_w' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_w' + filesuffix + ' width=300></td>\n')



    # Deal with r_prime
    # -------------------------------------------------------------------
    # For this we use the autoscale
    print '-------------------'
    field      = nc_file.variables['r_prime'][tstep,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    half_level = nc_file.variables['half_level'][:] / 1000.0
    print ' r_prime: ', field.shape
    #minlev = -0.04
    #maxlev =  0.04
    #con_levs = np.linspace(minlev, maxlev, 11)
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
    cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
    cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
    cax     = ax.contour (longs_v, half_level, field, colors='k')
    # Labels etc
    ax.set_title('r_prime\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_r_prime' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_r_prime' + filesuffix + ' width=300></td>\n')



    # Deal with b_prime
    # -------------------------------------------------------------------
    # For this we use the autoscale
    print '-------------------'
    field      = nc_file.variables['b_prime'][tstep,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    full_level = nc_file.variables['full_level'][:] / 1000.0
    print ' b_prime: ', field.shape
    #minlev = -0.04
    #maxlev =  0.04
    #con_levs = np.linspace(minlev, maxlev, 11)
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
    cax     = ax.contourf(longs_v, full_level, field, cmap=cmap)
    cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
    cax     = ax.contour (longs_v, full_level, field, colors='k')
    # Labels etc
    ax.set_title('b_prime\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_b_prime' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_b_prime' + filesuffix + ' width=300></td>\n')



    # Deal with the tracer
    # -------------------------------------------------------------------
    # For this we use the autoscale
    print '-------------------'
    field      = nc_file.variables['tracer'][tstep,:,:]
    longs_v    = nc_file.variables['longs_v'][:] / 1000.0
    half_level = nc_file.variables['half_level'][:] / 1000.0
    print ' tracer: ', field.shape
    #minlev = -0.04
    #maxlev =  0.04
    #con_levs = np.linspace(minlev, maxlev, 11)
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
    cmap    = cm.get_cmap('Greens', 11)
    cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
    cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
    ###cax     = ax.contour (longs_v, half_level, field, colors='k')
    # Labels etc
    ax.set_title('tracer\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
    ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
    ax.set_ylabel('Vertical distance (km)', fontsize=16)
    #plt.show()
    plt.savefig(plotdir + '/Plots/' + FileNamePart + '_tracer' + filesuffix, bbox_inches='tight')
    plt.close('all')

    if (output_type == 'web'):
      html_file.write ('<td><img src=Plots/' + FileNamePart + '_tracer' + filesuffix + ' width=300></td>\n')
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
