#!/usr/bin/env python
# -------------------------------------------------------------------
# Python code to read and plot the following quantities
#  vertical wind speed
#  effective buoyancy
#  geostrophic imbalance
#  hydrostatic imbalance
#  Source of vertical momentum flux
#  Horizontal divergence
#  Horizontal vorticity	
#  Zonal wind
#  Medidional wind
#  density pert (r_prime)
#  b_prime
#  tracer
#
# The plotting code is repeated for data in specified data directories
#
# Please edit the input details (e.g. location of data files)
# at the start of the main part of the code.
# This is located after the function definitions below.
#
# Ross Bannister, Dec 2016
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

# Set directories of the input data (each is a different experiment)
data_dir = []
data_dir.append('/home/data')
# Can add more directories in the same way

datafile = 'NewTruth.nc'


# As the code stands, the graphics are output to the same directories as the
# input data.  The graphics format is encapsulated postscript (eps).

# The time step of interest (0 for first timestep or in files of 3D data)
tstep    = 0

# Set the domain dimensions
nlongs   = 360
nlevs    = 60

# Output type - not yet implemented - will always output to web
# As there could be a large number of plots, there is an option to output results on a web page
# Choose 'eps' to output a collection of eps files (not on web page)
# Choose 'web' to output a collection of png images and an html script to view them on a web page
output_type = 'web'

if (output_type == 'web'):
  filesuffix = '.png'
else:
  filesuffix = '.eps'


for datadir in data_dir:

  plotdir = datadir
  print 'INPUT DIRECTORY  : ', datadir
  print 'PLOT DIRECTORY   : ', plotdir

  os.system('mkdir -p ' + plotdir + '/Plots')

  if (output_type == 'web'):
    # Set-up the html file
    html_file = open (plotdir + '/Plots.html', 'w')
    html_file.write ('<html>\n')
    html_file.write ('<h1>Model run plots</h1>\n')
    html_file.write (datadir)
    html_file.write ('<table>\n')
    html_file.write ('<tr>\n')


  # Deal with the zonal wind
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['u'][tstep,:,:]
  longs_u    = nc_file.variables['longs_u'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
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
#  cmap    = cm.get_cmap('Greys', 11)
  cmap    = cm.get_cmap('seismic', 11)
  cax     = ax.contourf(longs_u, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_u, half_level, field, colors='k')
  # Labels etc
  ax.set_title('u for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/u' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the meridional wind
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['v'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
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
#  cmap    = cm.get_cmap('Greys', 11)
  cmap    = cm.get_cmap('seismic', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  # Labels etc
  ax.set_title('v for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/v' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the vertical wind
  # -------------------------------------------------------------------
  # For this, we plot -0.4 to 0.4, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['w'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
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
  ax.set_title('w for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/w' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with r_prime
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['r_prime'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
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
#  cmap    = cm.get_cmap('Greys', 11)
  cmap    = cm.get_cmap('seismic', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  # Labels etc
  ax.set_title('rprime for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/r_prime' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')


  # Deal with b_prime
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['b_prime'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
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
#  cmap    = cm.get_cmap('Greys', 11)
  cmap    = cm.get_cmap('seismic', 11)
  cax     = ax.contourf(longs_v, full_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k')
  # Labels etc
  ax.set_title('bprime for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/b_prime' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')


  # Deal with the tracer
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['tracer'][tstep,:,:]
  field0     = nc_file.variables['tracer'][0,:,:]          # Initial tracer
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
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
#  cmap    = cm.get_cmap('Greys', 11)
  cmap    = cm.get_cmap('Greens', 11)
  cax     = ax.contourf(longs_v, half_level, field, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, colors='k')
  cax     = ax.contour (longs_v, half_level, field0, colors='k')
  # Labels etc
  ax.set_title('tracer for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/tracer' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the effective buoyancy
  # -------------------------------------------------------------------
  # For this we use the autoscale

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['b_effective'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' b_effective: ', field.shape
  minlev = -0.0005  # SET THESE ONLY FOR B+
  maxlev =  0.0015  # SET THESE ONLY FOR B+
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  #con_levs = np.linspace(minfield, maxfield, 9)
  con_levs = np.linspace(minlev, maxlev, 9)
  print con_levs
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
  cax     = ax.contourf(longs_v, full_level, field, levels=con_levs, cmap=cmap)
  cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k')
  # Labels etc
  ax.set_title('dbprime/dz + A^2 for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/b_eff' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the geostrophic imbalance
  # -------------------------------------------------------------------
  # For this we plot -1 to 1, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['geo_imbal'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' geo_imbal: ', field.shape
  minlev = -1.0
  maxlev =  1.0
  con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)

  matplotlib.rc('xtick', labelsize=16)
  matplotlib.rc('ytick', labelsize=16)
  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', len(con_levs))
  cax     = ax.contourf(longs_v, half_level, field, clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, half_level, field, clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('Geostrophic imbalance for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/geo_imbal' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the hydrostatic imbalance
  # -------------------------------------------------------------------
  # For this we plot -1 to 1, but mask out anything above 1, below -1

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['hydro_imbal'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' hydro_imbal: ', field.shape
  minlev = -1.0
  maxlev =  1.0
  con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  print 'Contour levels: ', con_levs
  # These two lines are tricks for making sure that the scales remain the same for all plots
  field[0][0] = minlev ; field[0][1] = maxlev
  field = remove_extremes (field, minlev, maxlev)

  matplotlib.rc('xtick', labelsize=16)
  matplotlib.rc('ytick', labelsize=16)
  fig, ax = plt.subplots()
  cmap    = cm.get_cmap('seismic', 9) #len(con_levs))
  cax     = ax.contourf(longs_v, full_level, field, cmap=cmap) #clevs=con_levs, cmap=cmap, extend='both')
  cbar    = fig.colorbar(cax, orientation='vertical', cmap=cmap) #ticks=con_levs, cmap=cmap)
  cax     = ax.contour (longs_v, full_level, field, colors='k') #clevs=con_levs, colors='k')
  # Labels etc
  ax.set_title('Hydrostatic imbalance for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/hydro_imbal' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')

  
  # Deal with the source of vertical momentum flux
  # -------------------------------------------------------------------
  # For this we use the autoscale, but with same size of +/- limits

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['wmom_source'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  print ' wmom_source: ', field.shape
  #minlev = -1.0
  #maxlev =  1.0
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  minlev   = -abs(max(abs(minfield), abs(maxfield)))
  maxlev   = abs(max(abs(minfield), abs(maxfield)))
  con_levs = np.linspace(minlev, maxlev, 11)
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs
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
  ax.set_title('Vertical momentum source for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/wmom_source' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the horizontal divergence
  # -------------------------------------------------------------------
  # For this we use the autoscale, but with same size of +/- limits

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['horiz_div'][tstep,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' horiz_div: ', field.shape
  #minlev = -1.0
  #maxlev =  1.0
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  minlev   = -abs(max(abs(minfield), abs(maxfield)))
  maxlev   = abs(max(abs(minfield), abs(maxfield)))
  con_levs = np.linspace(minlev, maxlev, 11)
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs
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
  ax.set_title('Horizontal divergence for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/horiz_div' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')



  # Deal with the horizontal vorticity
  # -------------------------------------------------------------------
  # For this we use the autoscale, but with same size of +/- limits

  nc_file = Dataset(datadir + '/' + datafile)
  print '-------------------'
  #print nc_file
  field      = nc_file.variables['horiz_vort'][tstep,:,:]
  longs_v    = nc_file.variables['longs_u'][:] / 1000.0
  full_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  print ' horiz_vort: ', field.shape
  #minlev = -1.0
  #maxlev =  1.0
  #con_levs = np.linspace(minlev, maxlev, 11)
  minfield = min_field(field)
  maxfield = max_field(field)
  print 'Min value of this field ', minfield
  print 'Max value of this field ', maxfield
  minlev   = -abs(max(abs(minfield), abs(maxfield)))
  maxlev   = abs(max(abs(minfield), abs(maxfield)))
  con_levs = np.linspace(minlev, maxlev, 11)
  mn       = np.mean(field)
  err      = field - mn
  rms      = np.sqrt(np.mean(err * err))
  print 'mean of this field      ', mn
  print 'RMS of this field       ', rms
  #print 'Contour levels: ', con_levs
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
  ax.set_title('Horizontal vorticity for t = ' + str(tstep*1800/3600) + 'h\nmin=' + str.format('%.2f' % minfield) + ', max=' + str.format('%.2f' % maxfield) + ', mean=' + str.format('%.5f' % mn) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
  ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
  ax.set_ylabel('Vertical distance (km)', fontsize=16)
  #plt.show()
  graphics_file_name = 'Plots/horiz_vort' + str(tstep) + filesuffix
  plt.savefig(plotdir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td><img src=' + graphics_file_name + ' width=300></td>\n')
