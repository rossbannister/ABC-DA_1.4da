import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from matplotlib import colors, cm
import matplotlib
import math


# ===================================================================
# ===================================================================
def max_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  maxlev  = []
  for lev in field:
    maxlev.append(max(lev))
  globalmax = max(maxlev)
  return globalmax

# ===================================================================
# ===================================================================
def min_field2d (field):
  # Subroutine to return maxmimum absolute value of field
  minlev  = []
  for lev in field:
    minlev.append(min(lev))
  globalmin = min(minlev)
  return globalmin

# ===================================================================
# ===================================================================
def max_field3d (field):
  # Subroutine to return maxmimum value of field
  maxt  = []
  for t in field:
    maxlev  = []
    for lev in t:
      maxlev.append(max(lev))
    maxt.append(max(maxlev))
  globalmax = max(maxt)
  return globalmax

# ===================================================================
# ===================================================================
def min_field3d (field):
  # Subroutine to return minimum value of field
  mint  = []
  for t in field:
    minlev  = []
    for lev in t:
      minlev.append(min(lev))
    mint.append(min(minlev))
  globalmin = min(mint)
  return globalmin

# ===================================================================
# ===================================================================
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
# ===================================================================
def flexi_round (value):
  s      = math.log10(value)
  si     = int(s-1.0)
  sm     = 10.0 ** float(si)
  modval = int(value/sm + 0.999999999999) * sm
  return modval


# ===================================================================
# ===================================================================
def OuterLoopBounds (J, Jb, Jo, outerloop, totalit):
  # To return an array of lines to indicate outer loop boundaries (for plots)
  minval3 = [np.nanmin(J[:-1]), np.nanmin(Jb[:-1]), np.nanmin(Jo[:-1])]
  minval  = np.nanmin(minval3)
  maxval3 = [np.nanmax(J[:-1]), np.nanmax(Jb[:-1]), np.nanmax(Jo[:-1])]
  maxval  = np.nanmax(maxval3)

  # Whenever the outer loop changes, create a line at the total iteration index
  loopboundsx = [[1.0, 1.0]]
  loopboundsy = [[minval, maxval]]
  for p0 in range(len(outerloop)-1):
    p1 = p0 + 1
    if (outerloop[p1] != outerloop[p0]):
      # x positions of the new line
      linex = [totalit[p1], totalit[p1]]
      # y positions of the new line
      liney = [minval, maxval]
      loopboundsx.append(linex)
      loopboundsy.append(liney)
  return loopboundsx, loopboundsy


# ===================================================================
# ===================================================================
def filterfield (field_in, lengths, scale):
  # Subroutine to return filtered field (bandpass down to the given scale)
  import numpy as np

  if scale == 1.0:
    # No filtering required
    print 'No filtering required'
    field_out = field_in
  else:
    # FFT the field
    field_fft = np.fft.fft(field_in)

    # Modify the field
    print 'Filtering the field'
    for k in range(0, len(lengths)):
      modfactor       = np.arctan(lengths[k] - scale) / np.pi + 0.5
      field_fft[:,k] *= modfactor

    # Inverse FFT the field
    field_out = np.real(np.fft.ifft(field_fft))
  return field_out


# ===================================================================
# ===================================================================
def calc_geoimbal (rp, v, dist, f, C, nlongs, nlevs):
  # Calculate the geostrophic imbalance diagnostic
  import numpy as np

  # Set-up
  recipdx = 1.0 / (dist[1] - dist[0])
  gi_1     = np.zeros((nlevs,nlongs))
  gi_2     = np.zeros((nlevs,nlongs))

  # Calculate unnormalised geostrophic imbalance, terms 1 and 2
  for z in range(0,nlevs):
    for x in range(1,nlongs):
      gi_1[z,x] = C * recipdx * (rp[z,x] - rp[z,x-1])
      gi_2[z,x] = -1.0 * f * (v[z,x] + v[z,x-1]) / 2.0

  # Normalise
  # Calculate the rms of term 1
  mean    = np.mean(gi_1[:,1:nlongs])
  rms_1   = np.sqrt(np.mean( (gi_1[:,1:nlongs]-mean)*(gi_1[:,1:nlongs]-mean) ))
  # Calculate the rms of term 2
  mean    = np.mean(gi_2[:,1:nlongs])
  rms_2   = np.sqrt(np.mean( (gi_2[:,1:nlongs]-mean)*(gi_2[:,1:nlongs]-mean) ))
  gi      = (gi_1 + gi_2) / (rms_1 + rms_2)

  return gi


# ===================================================================
# ===================================================================
def calc_hydimbal (rp, bp, full_lev, half_lev, C, nlongs, nlevs):
  # Calculate the hydrostatic imbalance diagnostic
  import numpy as np

  # Set-up
  hi_1 = np.zeros((nlevs,nlongs))
  hi_2 = np.zeros((nlevs,nlongs))

  # Calculate unnormalised hydrostatic imbalance, terms 1 and 2
  for x in range(0,nlongs):
    for z in range(0,nlevs-1):
      hi_1[z,x] = C * (rp[z+1,x] - rp[z,x]) / (half_lev[z+1] - half_lev[z])
      hi_2[z,x] = -1.0 * bp[z,x]

  # Normalise
  # Calculate the rms of term 1
  mean    = np.mean(hi_1[0:nlevs,:])
  rms_1   = np.sqrt(np.mean( (hi_1[0:nlevs,:]-mean)*(hi_1[0:nlevs,:]-mean) ))
  # Calculate the rms of term 2
  mean    = np.mean(hi_2[0:nlevs,:])
  rms_2   = np.sqrt(np.mean( (hi_2[0:nlevs,:]-mean)*(hi_2[0:nlevs,:]-mean) ))
  hi      = (hi_1 + hi_2) / (rms_1 + rms_2)

  return hi


# ===================================================================
# ===================================================================
def balance_scale (bal_file1, bal_file2, t, hori_lens, scalelim, C_param, f, nlongs, nlevs):
  # 1. Read-in fields from file(s)
  # 2. Find their difference (if appropriate)
  # 3. Degrade their resolution
  # 4. Compute imbalance diagnostics
  # NOTE: If bal_file2 is an empty name then we plot data from bal_file1
  #       If bal_file2 is not empty then we plot data from bal_file2 - bal_file1

  # Read-in the relevant fields from file 1
  # ---------------------------------------
  nc_file    = Dataset(bal_file1)
  bp1        = nc_file.variables['b_prime'][t,:,:]
  rp1        = nc_file.variables['r_prime'][t,:,:]
  v1         = nc_file.variables['v'][t,:,:]
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close

  ##### Compute the balances associated with file 1
  gimbal_rms_scales1 = np.zeros(len(scalelim))
  himbal_rms_scales1 = np.zeros(len(scalelim))
  for scale in range (0,len(scalelim)):
    print '    Dealing with lengthscale band (file 1) ', scalelim[scale]
    ### Filter the fields to remove specific scales
    print '    Filtering fields'
    bp_filt = filterfield(bp1, hori_lens, scalelim[scale])
    rp_filt = filterfield(rp1, hori_lens, scalelim[scale])
    v_filt  = filterfield(v1,  hori_lens, scalelim[scale])
    ### Calculate the diagnostic fields
    gimbal  = calc_geoimbal(rp_filt, v_filt, 1000.0*longs_v, f, C_param, nlongs, nlevs)
    himbal  = calc_hydimbal(rp_filt, bp_filt, 1000.0*full_level, 1000.0*half_level, C_param, nlongs, nlevs)
    ### Compute the rms values
    gimbal_rms = np.sqrt(np.mean( gimbal[:,1:nlongs]*gimbal[:,1:nlongs]) )
    himbal_rms = np.sqrt(np.mean( himbal[0:nlevs,:]*himbal[0:nlevs,:] ))
    ### Store the values
    gimbal_rms_scales1[scale] = gimbal_rms
    himbal_rms_scales1[scale] = himbal_rms


  # Set-up data structures even if the contents are not calculated
  # --------------------------------------------------------------
  gimbal_rms_scales2     = np.zeros(len(scalelim))
  himbal_rms_scales2     = np.zeros(len(scalelim))
  gimbal_rms_scales_diff = np.zeros(len(scalelim))
  himbal_rms_scales_diff = np.zeros(len(scalelim))


  # Read-in the relevant fields from file 2 (if appropriate)
  # --------------------------------------------------------
  if len(bal_file2) > 0:
    nc_file   = Dataset(bal_file2)
    bp2       = nc_file.variables['b_prime'][t,:,:]
    rp2       = nc_file.variables['r_prime'][t,:,:]
    v2        = nc_file.variables['v'][t,:,:]
    nc_file.close

    ##### Compute the balances associated with file 2
    for scale in range (0,len(scalelim)):
      print '    Dealing with lengthscale band (file 2) ', scalelim[scale]
      ### Filter the fields to remove specific scales
      print '    Filtering fields'
      bp_filt = filterfield(bp2, hori_lens, scalelim[scale])
      rp_filt = filterfield(rp2, hori_lens, scalelim[scale])
      v_filt  = filterfield(v2,  hori_lens, scalelim[scale])
      ### Calculate the diagnostic fields
      gimbal  = calc_geoimbal(rp_filt, v_filt, 1000.0*longs_v, f, C_param, nlongs, nlevs)
      himbal  = calc_hydimbal(rp_filt, bp_filt, 1000.0*full_level, 1000.0*half_level, C_param, nlongs, nlevs)
      ### Compute the rms values
      gimbal_rms = np.sqrt(np.mean( gimbal[:,1:nlongs]*gimbal[:,1:nlongs]) )
      himbal_rms = np.sqrt(np.mean( himbal[0:nlevs,:]*himbal[0:nlevs,:] ))
      ### Store the values
      gimbal_rms_scales2[scale] = gimbal_rms
      himbal_rms_scales2[scale] = himbal_rms


    # Deal with the differences between the two sets of fields (file 2 - file 1)
    # Take the difference
    bp_diff = bp2 - bp1
    rp_diff = rp2 - rp1
    v_diff  = v2 - v1

    for scale in range (0,len(scalelim)):
      print '    Dealing with lengthscale band (difference) ', scalelim[scale]
      ### Filter the fields to remove specific scales
      print '    Filtering fields'
      bp_filt = filterfield(bp_diff, hori_lens, scalelim[scale])
      rp_filt = filterfield(rp_diff, hori_lens, scalelim[scale])
      v_filt  = filterfield(v_diff,  hori_lens, scalelim[scale])
      ### Calculate the diagnostic fields
      gimbal  = calc_geoimbal(rp_filt, v_filt, 1000.0*longs_v, f, C_param, nlongs, nlevs)
      himbal  = calc_hydimbal(rp_filt, bp_filt, 1000.0*full_level, 1000.0*half_level, C_param, nlongs, nlevs)
      ### Compute the rms values
      gimbal_rms = np.sqrt(np.mean( gimbal[:,1:nlongs]*gimbal[:,1:nlongs]) )
      himbal_rms = np.sqrt(np.mean( himbal[0:nlevs,:]*himbal[0:nlevs,:] ))
      ### Store the values
      gimbal_rms_scales_diff[scale] = gimbal_rms
      himbal_rms_scales_diff[scale] = himbal_rms

  return gimbal_rms_scales1, himbal_rms_scales1, gimbal_rms_scales2, himbal_rms_scales2, gimbal_rms_scales_diff, himbal_rms_scales_diff




# ===================================================================
# ===================================================================
def ReadObs_file (filename):
  # Read-in observations file

  u_obs         = []   ;   u_modelobs         = []   ;   u_trueobs         = []
  v_obs         = []   ;   v_modelobs         = []   ;   v_trueobs         = []
  w_obs         = []   ;   w_modelobs         = []   ;   w_trueobs         = []
  r_obs         = []   ;   r_modelobs         = []   ;   r_trueobs         = []
  b_obs         = []   ;   b_modelobs         = []   ;   b_trueobs         = []
  tracer_obs    = []   ;   tracer_modelobs    = []   ;   tracer_trueobs    = []
  horizwind_obs = []   ;   horizwind_modelobs = []   ;   horizwind_trueobs = []
  totalwind_obs = []   ;   totalwind_modelobs = []   ;   totalwind_trueobs = []
  
  obsfile = open (filename, 'r')
  line    = obsfile.readline()
  line    = obsfile.readline()
  version = int(line[17:])
  for loop in range(7):
    line    = obsfile.readline()

  while (line != ''):
    for loop in range(11):
      line    = obsfile.readline()
    obs_of = int(line[17:])
    for loop in range(2):
      line    = obsfile.readline()
    truth_available = line[18:19]
    line     = obsfile.readline()
    if truth_available == 'T':
      truth  = float(line[17:])
    else:
      truth  = 999.0
    line     = obsfile.readline()
    obs      = float(line[17:])
    line     = obsfile.readline()
    stddev   = float(line[17:])
    line     = obsfile.readline()
    modelobs = float(line[17:])
    for loop in range(5):
      line    = obsfile.readline()
    # print truth_available, truth, obs, modelobs
    if obs_of == 1:
      # Zonal wind
      u_obs.append(obs)
      u_modelobs.append(modelobs)
      u_trueobs.append(truth)
    elif obs_of ==2:
      # Meridional wind
      v_obs.append(obs)
      v_modelobs.append(modelobs)
      v_trueobs.append(truth)
    elif obs_of ==3:
      # Vertical wind
      w_obs.append(obs)
      w_modelobs.append(modelobs)
      w_trueobs.append(truth)
    elif obs_of ==4:
      # density pert
      r_obs.append(obs)
      r_modelobs.append(modelobs)
      r_trueobs.append(truth)
    elif obs_of ==5:
      # buoyancy
      b_obs.append(obs)
      b_modelobs.append(modelobs)
      b_trueobs.append(truth)
    elif obs_of ==6:
      # tracer
      tracer_obs.append(obs)
      tracer_modelobs.append(modelobs)
      tracer_trueobs.append(truth)
    elif obs_of ==7:
      # horizontal wind
      horizwind_obs.append(obs)
      horizwind_modelobs.append(modelobs)
      horizwind_trueobs.append(truth)
    elif obs_of ==8:
      # total wind
      totalwind_obs.append(obs)
      totalwind_modelobs.append(modelobs)
      totalwind_trueobs.append(truth)

  obsfile.close

  return u_obs, v_obs, w_obs, r_obs, b_obs, tracer_obs, horizwind_obs, totalwind_obs, u_modelobs, v_modelobs, w_modelobs, r_modelobs, b_modelobs, tracer_modelobs, horizwind_modelobs, totalwind_modelobs, u_trueobs, v_trueobs, w_trueobs, r_trueobs, b_trueobs, tracer_trueobs, horizwind_trueobs, totalwind_trueobs
  
  
  
# ===================================================================
# ===================================================================
def plotfields (filename, prefix, output_type, html_file, plotdir, TimeOutput):

# TimeOutput:
# Plot all times "AllTimes" or just first and last times "FirstLast" in above options?

  if (output_type == 'web'):
    filesuffix = '.png'
  else:
    filesuffix = '.eps'


  # Deal with the zonal wind
  # -------------------------------------------------------------------
  nc_file = Dataset(filename)
  field   = nc_file.variables['u'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<table>')
    html_file.write ('<tr><td>u</td></tr>\n')
    html_file.write ('<tr>\n')

  # Find how many time steps are in this file
  shape      = field.shape
  ntimes     = shape[0]
  
  longs_u    = nc_file.variables['longs_u'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  times      = nc_file.variables['time'][:] / 3600.0
  nc_file.close
  #minlev = -0.4
  #maxlev =  0.4
  minfield = min_field3d(field)
  maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)
  
  # ----- Determine the times that we should print at (need to do this only for u field)
  if TimeOutput == 'AllTimes':
    # Plot for all times
    PlotTimeList = range(ntimes)
  elif TimeOutput == 'FirstLast':
    PlotTimeList = [0,ntimes-1]
  else:
    PlotTimeList = []
  # -----

  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting u at time', t
      fieldt = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_u, half_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_u, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('u for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_u_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')


  # Deal with the meridional wind
  # -------------------------------------------------------------------
  nc_file = Dataset(filename)
  field   = nc_file.variables['v'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>v</td></tr>\n')
    html_file.write ('<tr>\n')

  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  #minlev = -0.4
  #maxlev =  0.4
  minfield = min_field3d(field)
  maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting v at time', t
      fieldt = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, half_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('v for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_v_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')


  # Deal with the vertical wind
  # -------------------------------------------------------------------
  # For this, we plot -0.4 to 0.4, but mask out anything above 1, below -1
  nc_file = Dataset(filename)
  field   = nc_file.variables['w'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>w</td></tr>\n')
    html_file.write ('<tr>\n')

  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  minfield = -0.4
  maxfield =  0.4
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting w at time', t
      fieldt   = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        # These two lines are tricks for making sure that the scales remain the same for all plots
        fieldt[0][0] = minlev ; fieldt[0][1] = maxlev
        fieldt = remove_extremes (fieldt, minfield, maxfield)
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, full_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, full_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('w for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_w_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')




  # Deal with the r_prime
  # -------------------------------------------------------------------
  nc_file = Dataset(filename)
  field   = nc_file.variables['r_prime'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>r_prime</td></tr>\n')
    html_file.write ('<tr>\n')

  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  #minlev = -0.4
  #maxlev =  0.4
  minfield = min_field3d(field)
  maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting r_prime at time', t
      fieldt = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, half_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('r_prime for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_rp_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')




  # Deal with b_prime
  # -------------------------------------------------------------------
  nc_file = Dataset(filename)
  field   = nc_file.variables['b_prime'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>b_prime</td></tr>\n')
    html_file.write ('<tr>\n')

  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  full_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  #minlev = -0.4
  #maxlev =  0.4
  minfield = min_field3d(field)
  maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting b_prime at time', t
      fieldt = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, full_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('b_prime for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_bp_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')



  # Deal with the tracer
  # -------------------------------------------------------------------
  nc_file = Dataset(filename)
  field   = nc_file.variables['tracer'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>tracer</td></tr>\n')
    html_file.write ('<tr>\n')

  # Find how many time steps are in this file
  shape      = field.shape
  ntimes     = shape[0]
  
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  #minlev = -0.4
  #maxlev =  0.4
  minfield = min_field3d(field)
  maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting tracer at time', t
      fieldt = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('Greens', len(con_levs))
        cax     = ax.contourf(longs_v, half_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('tracer for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_tracer_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')



  # Deal with the geostrophic balance
  # -------------------------------------------------------------------
  # For this, we plot -1.0 to 1.0
  nc_file = Dataset(filename)
  field   = nc_file.variables['geo_imbal'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>geo imbal</td></tr>\n')
    html_file.write ('<tr>\n')
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['half_level'][:] / 1000.0
  nc_file.close
  minfield = -1.0
  maxfield =  1.0
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting geo_imbal at time', t
      fieldt   = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        # These two lines are tricks for making sure that the scales remain the same for all plots
        fieldt[0][0] = minlev ; fieldt[0][1] = maxlev
        fieldt = remove_extremes (fieldt, minfield, maxfield)
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, half_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, half_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('geo_imbal for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_geoim_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')



  # Deal with the hydrostatic balance
  # -------------------------------------------------------------------
  # For this, we plot -1.0 to 1.0
  nc_file = Dataset(filename)
  field   = nc_file.variables['hydro_imbal'][:,:,:]

  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<tr><td>hydro_imbal</td></tr>\n')
    html_file.write ('<tr>\n')
  longs_v    = nc_file.variables['longs_v'][:] / 1000.0
  half_level = nc_file.variables['full_level'][:] / 1000.0
  nc_file.close
  minfield = -1.0
  maxfield =  1.0
  con_levs = np.linspace(minfield, maxfield, 11)
  for t in range(ntimes):
    if (t in PlotTimeList):
      print 'Plotting hydro_imbal at time', t
      fieldt   = field[t][:][:]
      minlev = min_field2d(fieldt)
      maxlev = max_field2d(fieldt)
      rms    = np.sqrt(np.mean(fieldt * fieldt))
      print 'Min value of this field ', minlev
      print 'Max value of this field ', maxlev
      print 'RMS of this field       ', rms
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        print 'Contour levels: ', con_levs
        # These two lines are tricks for making sure that the scales remain the same for all plots
        fieldt[0][0] = minlev ; fieldt[0][1] = maxlev
        fieldt = remove_extremes (fieldt, minfield, maxfield)
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs_v, full_level, fieldt, clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs_v, full_level, fieldt, clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title('hydro_imbal for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % rms), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_hydroim_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>\n')
    html_file.write ('</table>')

  return



# ===================================================================
# ===================================================================
def plot_1_field_diff (filename1, filename2, variable, x_grid_name, z_grid_name, prefix, filesuffix, output_type, html_file, plotdir, TimeOutput):
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables[variable][:,:,:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables[variable][:,:,:]
  field    = field2[:,:,:] - field1[:,:,:]
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('<table>')
    html_file.write ('<tr><td>' + variable + ' difference</td></tr>\n')
    html_file.write ('<tr>\n')
  # Find how many time steps are in this file
  shape    = field.shape
  ntimes   = shape[0]
  longs    = nc_file1.variables[x_grid_name][:] / 1000.0
  levs     = nc_file1.variables[z_grid_name][:] / 1000.0
  times    = nc_file1.variables['time'][:] / 3600.0
  nc_file1.close
  nc_file2.close
  # Minimum and maximum levels of the difference field
  if (variable == 'geo_imbal' or variable == 'hydro_imbal'):
    minfield = -1.0
    maxfield = 1.0
  else:
    minfield = min_field3d(field)
    maxfield = max_field3d(field)
  con_levs = np.linspace(minfield, maxfield, 11)

  # ----- Determine the times that we should print at
  if TimeOutput == 'AllTimes':
    # Plot for all times
    PlotTimeList = range(ntimes)
  elif TimeOutput == 'FirstLast':
    PlotTimeList = [0,ntimes-1]
  else:
    PlotTimeList = []
  # -----

  RMS1      = []
  Mean1     = []
  RMS2      = []
  Mean2     = []
  RMSd      = []
  Meand     = []
  for t in range(ntimes):
    print 'Dealing with u at time', t
    # Various averages of each field separately
    RMS1.append (np.sqrt(np.mean(field1[t][:][:] * field1[t][:][:])))
    Mean1.append(np.mean(field1[t][:][:]))
    RMS2.append (np.sqrt(np.mean(field2[t][:][:] * field2[t][:][:])))
    Mean2.append(np.mean(field2[t][:][:]))
    # Deal with the differences, field2 - field1
    minlev = min_field2d(field[t][:][:])
    maxlev = max_field2d(field[t][:][:])
    print 'Min value of this field ', minlev
    print 'Max value of this field ', maxlev
    # Compute the RMS and mean of the difference
    RMSd.append (np.sqrt(np.mean(field[t][:][:] * field[t][:][:])))
    Meand.append(np.mean(field[t][:][:]))
    
    if (t in PlotTimeList):
      if (minlev == maxlev):
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td>Constant value ' + str(minlev) + '</td>\n')
      else:
        #print 'Contour levels: ', con_levs
        matplotlib.rc('xtick', labelsize=16)
        matplotlib.rc('ytick', labelsize=16)
        fig, ax = plt.subplots()
        cmap    = cm.get_cmap('seismic', len(con_levs))
        cax     = ax.contourf(longs, levs, field[t][:][:], clevs=con_levs, cmap=cmap, extend='both')
        cbar    = fig.colorbar(cax, orientation='vertical', ticks=con_levs, cmap=cmap)
        cax     = ax.contour (longs, levs, field[t][:][:], clevs=con_levs, colors='k')
        # Labels etc
        ax.set_title(variable + ' for t = ' + str(times[t]) + 'h\nmin=' + str.format('%.2f' % minlev) + ', max=' + str.format('%.2f' % maxlev) + ', rms=' + str.format('%.5f' % RMSd[t]), fontsize=16)
        ax.set_xlabel('Longitudinal distance (km)', fontsize=16)
        ax.set_ylabel('Vertical distance (km)', fontsize=16)
        #plt.show()
        figfilename = 'Plots/' + prefix + '_' + variable + '_' + str(t) + filesuffix
        plt.savefig(plotdir + '/' + figfilename, bbox_inches='tight')
        if (output_type == 'web' and TimeOutput != 'None'):
          html_file.write ('<td><img src=' + figfilename  + ' width=300></td>\n')
        plt.close('all')
  # Store the RMS time sequence for this variable
  if (output_type == 'web' and TimeOutput != 'None'):
    html_file.write ('</tr>')
    html_file.write ('</table>')
  return RMS1, Mean1, RMS2, Mean2, RMSd, Meand



# ===================================================================
# ===================================================================
def plot_diff_fields (filename1, filename2, prefix, output_type, html_file, plotdir, C_param, f, TimeOutput):

  # To plot difference fields
  # Also to return:
  #   Mean of data in filename1 (for each field)
  #   RMS of data in filename1 (for each field)
  #   Mean of data in filename2 (for each field)
  #   RMS of data in filename2 (for each field)
  #   Mean of difference of data in filename2 - filename1 (for each field)
  #   RMS of difference of data in filename2 - filename1 (for each field)
  # Also to return the difference of energies
  # Also to return the difference of scale-dependent imbalances

  # TimeOutput:
  # Plot all times "AllTimes" or just first and last times "FirstLast" in above options?
  # The option "None" is also used to turn-off plotting altogether (and just collect data).


  # Set-up the structures containing the mean and RMS stuff
  # This will be a list of sub-lists
  # Each sub-list corresponds to a time sequence of RMS/mean values/differences
  # The main list is for each quantity in order processed
  Mean_file1    = []
  RMS_file1     = []
  Mean_file2    = []
  RMS_file2     = []
  Mean_diff     = []
  RMS_diff      = []
  QuantityNames = []

  if (output_type == 'web'):
    filesuffix = '.png'
  else:
    filesuffix = '.eps'


  # Deal with the zonal wind
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'u', 'longs_u', 'half_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('u')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the meridional wind
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'v', 'longs_v', 'half_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('v')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the vertical wind
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'w', 'longs_v', 'full_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('w')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the r_prime
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'r_prime', 'longs_v', 'half_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('r_prime')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with b_prime
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'b_prime', 'longs_v', 'full_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('b_prime')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the tracer
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'tracer', 'longs_v', 'half_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('tracer')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the geostrophic imbalance
  # -------------------------------------------------------------------
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'geo_imbal', 'longs_v', 'half_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('geo_imbal')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)


  # Deal with the hydrostatic imbalance
  # -------------------------------------------------------------------
  # For this, we plot -1.0 to 1.0
  RMS1, Mean1, RMS2, Mean2, RMSd, Meand = plot_1_field_diff (filename1, filename2, 'hydro_imbal', 'longs_v', 'full_level', prefix, filesuffix, output_type, html_file, plotdir, TimeOutput)
  # Append the time sequence to the lists that will contain all quantities
  QuantityNames.append('hydro_imbal')
  Mean_file1.append(Mean1)
  RMS_file1.append(RMS1)
  Mean_file2.append(Mean2)
  RMS_file2.append(RMS2)
  Mean_diff.append(Meand)
  RMS_diff.append(RMSd)
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables['hydro_imbal'][:,:,:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables['hydro_imbal'][:,:,:]
  field    = field2[:,:,:] - field1[:,:,:]



  # Also obtain (but do not plot) energy data
  # -------------------------------------------------------------------
  # Kinetic energy
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables['ke'][:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables['ke'][:]
  field    = field2[:] - field1[:]
  nc_file1.close
  nc_file2.close
  QuantityNames.append('ke')
  RMS_file1.append(np.abs(field1))
  Mean_file1.append(field1)
  RMS_file2.append(np.abs(field2))
  Mean_file2.append(field2)
  RMS_diff.append(np.abs(field))
  Mean_diff.append(field)

  # Buoyant energy
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables['be'][:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables['be'][:]
  field    = field2[:] - field1[:]
  nc_file1.close
  nc_file2.close
  QuantityNames.append('be')
  RMS_file1.append(np.abs(field1))
  Mean_file1.append(field1)
  RMS_file2.append(np.abs(field2))
  Mean_file2.append(field2)
  RMS_diff.append(np.abs(field))
  Mean_diff.append(field)

  # Elastic energy
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables['ee'][:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables['ee'][:]
  field    = field2[:] - field1[:]
  nc_file1.close
  nc_file2.close
  QuantityNames.append('ee')
  RMS_file1.append(np.abs(field1))
  Mean_file1.append(field1)
  RMS_file2.append(np.abs(field2))
  Mean_file2.append(field2)
  RMS_diff.append(np.abs(field))
  Mean_diff.append(field)

  # Total energy
  nc_file1 = Dataset(filename1)
  field1   = nc_file1.variables['te'][:]
  nc_file2 = Dataset(filename2)
  field2   = nc_file2.variables['te'][:]
  field    = field2[:] - field1[:]
  nc_file1.close
  nc_file2.close
  QuantityNames.append('te')
  RMS_file1.append(np.abs(field1))
  Mean_file1.append(field1)
  RMS_file2.append(np.abs(field2))
  Mean_file2.append(field2)
  RMS_diff.append(np.abs(field))
  Mean_diff.append(field)


  # Also obtain (but do not plot) scale-dependent imbalances
  # -------------------------------------------------------------------

  # Determine the number of longitudes and levels
  nc_file1 = Dataset(filename1)
  nlongs   = len(nc_file1.variables['longs_u'][:])
  nlevs    = len(nc_file1.variables['half_level'][:])
  times    = nc_file1.variables['time'][:]
  ntimes   = len(times)
  nc_file1.close
  Lx       = 1.5 * float(nlongs)

  # Wavenumbers and scales
  hori_wns      = np.linspace(0,nlongs-1,nlongs)
  hori_lens     = np.zeros(nlongs)
  hori_lens[1:nlongs/2+1] = Lx / hori_wns[1:nlongs/2+1]
  hori_lens[0]  = Lx
  for l in range(1, nlongs/2):
    hori_lens[nlongs-l] = hori_lens[l]

  # Define the scales of interest (km)
  scalelim = []
  scalelim.append(100.0)
  scalelim.append(10.0)
  scalelim.append(1.0)

  # Compute the scale-dependent balances
  RMS_gim1     = []   ; RMS_him1     = []
  RMS_gim2     = []   ; RMS_him2     = []
  RMS_gim_diff = []   ; RMS_him_diff = []
  for t in range(ntimes):
    gimbal_rms_scales1, himbal_rms_scales1, gimbal_rms_scales2, himbal_rms_scales2, gimbal_rms_scales_diff, himbal_rms_scales_diff = balance_scale (filename1, filename2, t, hori_lens, scalelim, C_param, f, nlongs, nlevs)
    RMS_gim1.append(gimbal_rms_scales1)
    RMS_him1.append(himbal_rms_scales1)
    RMS_gim2.append(gimbal_rms_scales2)
    RMS_him2.append(himbal_rms_scales2)
    RMS_gim_diff.append(gimbal_rms_scales_diff)
    RMS_him_diff.append(himbal_rms_scales_diff)

  # The above structures RMS_gim and RMS_him are time, scale index.  We want scale index, time.
  # Create new structures with the required order
  # Geostrophic balance at each scale limit

  for s in range(len(scalelim)):
    RMS_gim1_temp = []      ;   RMS_him1_temp = []
    RMS_gim2_temp = []      ;   RMS_him2_temp = []
    RMS_gim_diff_temp = []  ;   RMS_him_diff_temp = []
    for t in range(ntimes):
      RMS_gim1_temp.append(RMS_gim1[t][s])
      RMS_him1_temp.append(RMS_him1[t][s])
      RMS_gim2_temp.append(RMS_gim2[t][s])
      RMS_him2_temp.append(RMS_him2[t][s])
      RMS_gim_diff_temp.append(RMS_gim_diff[t][s])
      RMS_him_diff_temp.append(RMS_him_diff[t][s])
    QuantityNames.append('GeoImbal_' + str(scalelim[s]))
    QuantityNames.append('HydroImbal_' + str(scalelim[s]))
    RMS_file1.append(RMS_gim1_temp)
    RMS_file1.append(RMS_him1_temp)
    RMS_file2.append(RMS_gim2_temp)
    RMS_file2.append(RMS_him2_temp)
    RMS_diff.append(RMS_gim_diff_temp)
    RMS_diff.append(RMS_him_diff_temp)

  if (output_type == 'web'):
    html_file.write ('</tr>\n')

  return times, QuantityNames, Mean_file1, Mean_file2, Mean_diff, RMS_file1, RMS_file2, RMS_diff




# ===================================================================
# ===================================================================
def dump_scalar_time_seq (times, quantities, code, TRUTH_mean, DATA_mean, DATA_err, TRUTH_RMS, DATA_RMS, DATA_RMSE, output_dir):
  # Subroutine to dump time sequencies of mean_diff, and rms_diff data to a file
  # These structures are lists of sublists.
  # Each sublist is the time sequence for each quantity
  # The scale-resolved balance diagnostics have only rms_diff entries

  Num_rmss  = len(DATA_RMSE)
  Num_means = len(DATA_mean)

  output_file = open (output_dir + '/' + code + '.dat', 'w')
  output_file.write (code + ' values and error data\n')
  output_file.write ('===================================\n')
  output_file.write ('time\n')
  output_string = ''
  for number in times:
    output_string += str(number) + '  '
  output_file.write (output_string + '\n')

  # Go over each item that has a mean value
  for item in range(len(DATA_mean)):
    output_file.write (quantities[item] + '\n')
    # The mean value of the truth
    output_string = ''
    for number in TRUTH_mean[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The mean value of the data
    output_string = ''
    for number in DATA_mean[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The mean value of the error in the data
    output_string = ''
    for number in DATA_err[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The rms value of the truth
    output_string = ''
    for number in TRUTH_RMS[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The rms value of the data
    output_string = ''
    for number in DATA_RMS[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The rms value of the error in the data
    output_string = ''
    for number in DATA_RMSE[item]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')

  # Go over each extra item that has an rms value only
  for item in range(Num_rmss - Num_means):
    output_file.write (quantities[item+Num_means] + '\n')
    output_string = ''
    for number in TRUTH_RMS[item+Num_means]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The rms value of the data
    output_string = ''
    for number in DATA_RMS[item+Num_means]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
    # The rms value of the error in the data
    output_string = ''
    for number in DATA_RMSE[item+Num_means]:
      output_string += str(number) + '  '
    output_file.write (output_string + '\n')
  output_file.close()
  return


# ===================================================================
# ===================================================================
def plot_scalar_time_seq (time, cycle_bound_times, quantities, code, TRUTH_mean, DATA_mean, DATA_err, TRUTH_RMS, DATA_RMS, FREE_RMS, DATA_RMSE, FREE_RMSE, output_type, html_file, plot_dir):
  # Subroutine to plot time sequencies of mean_diff, and rms_diff
  # These structures are lists of sublists.
  # Each sublist is the time sequence for each quantity
  # The scale-resolved balance diagnostics have only rms_diff entries

  if (output_type == 'web'):
    filesuffix = '.png'
    html_file.write ('<table>\n')
  else:
    filesuffix = '.eps'


  # ===== Plot the truth =====
  # Plot truth only if code is 'bg'
  if code == 'bg':
    # All wind components, truth
    # --------------------------
    print 'True winds'
    if (output_type == 'web'):
      html_file.write ('<tr><td>Wind component truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m/s)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[0][:],  linewidth=1, ls='solid', color='red', label='u RMS')
    ax.plot(time[:], TRUTH_RMS[1][:],  linewidth=1, ls='solid', color='blue', label='v RMS')
    ax.plot(time[:], TRUTH_RMS[2][:],  linewidth=1, ls='solid', color='green', label='w RMS')
    ax.plot(time[:], TRUTH_mean[0][:], linewidth=2, ls='dotted', color='red', label='u mean')
    ax.plot(time[:], TRUTH_mean[1][:], linewidth=2, ls='dotted', color='blue', label='v mean')
    ax.plot(time[:], TRUTH_mean[2][:], linewidth=2, ls='dotted', color='green', label='w mean')
    ax.legend(loc='best') #(loc='upper right')
    plt.title('RMS and mean true values of wind components')
    graphics_file_name = 'Plots/Wind_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # rp, truth
    # -------------------
    print 'True rp'
    if (output_type == 'web'):
      html_file.write ('<tr><td>rp truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[3][:], linewidth=1, ls='solid', color='red', label='RMS')
    ax.plot(time[:], TRUTH_mean[3][:], linewidth=2, ls='dotted', color='red', label='mean')
    plt.title('RMS and mean true values of rp')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/rp_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # bp, truth
    # -------------------
    print 'True bp'
    if (output_type == 'web'):
      html_file.write ('<tr><td>bp truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m/s2)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[4][:], linewidth=1, ls='solid', color='blue', label='RMS')
    ax.plot(time[:], TRUTH_mean[4][:], linewidth=2, ls='dotted', color='blue', label='mean')
    plt.title('RMS and mean true values of buoyancy')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/bp_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # tracer, truth
    # -------------------
    print 'True tracer'
    if (output_type == 'web'):
      html_file.write ('<tr><td>tracer truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[5][:], linewidth=1, ls='solid', color='blue', label='RMS')
    ax.plot(time[:], TRUTH_mean[5][:], linewidth=2, ls='dotted', color='blue', label='mean')
    plt.title('RMS and mean true values of tracer')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/tracer_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # imbalances, truth
    # -------------------
    print 'True imbalances'
    if (output_type == 'web'):
      html_file.write ('<tr><td>imbalances truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[6][:], linewidth=1, ls='solid', color='blue', label='Geo imbal RMS')
    ax.plot(time[:], TRUTH_RMS[7][:], linewidth=1, ls='solid', color='orange', label='Hyd imbal RMS')
    ax.plot(time[:], TRUTH_mean[6][:], linewidth=2, ls='dotted', color='blue', label='Hyd imbal mean')
    ax.plot(time[:], TRUTH_mean[7][:], linewidth=2, ls='dotted', color='orange', label='Geo imbal mean')
    plt.title('RMS and mean true values of imbalances')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/imbal_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # energy, truth
    # -------------------
    print 'True energy'
    if (output_type == 'web'):
      html_file.write ('<tr><td>energy truth<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m2/s2)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], TRUTH_RMS[8][:], linewidth=2, ls='solid', color='red', label='KE RMS')
    ax.plot(time[:], TRUTH_RMS[9][:], linewidth=2, ls='solid', color='blue', label='BE RMS')
    ax.plot(time[:], TRUTH_RMS[10][:], linewidth=2, ls='solid', color='green', label='EE RMS')
    ax.plot(time[:], TRUTH_RMS[11][:], linewidth=2, ls='solid', color='black', label='TE RMS')
    ax.plot(time[:], TRUTH_mean[8][:], linewidth=2, ls='dotted', color='red', label='KE mean')
    ax.plot(time[:], TRUTH_mean[9][:], linewidth=2, ls='dotted', color='blue', label='BE mean')
    ax.plot(time[:], TRUTH_mean[10][:], linewidth=2, ls='dotted', color='green', label='EE mean')
    ax.plot(time[:], TRUTH_mean[11][:], linewidth=2, ls='dotted', color='black', label='TE mean')
    plt.title('RMS and mean true values of energy')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/energy_truth' + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # Scale-resolved balances, truth (RMS only)
    # -------------------
    print 'True scale resolved balances'
    if (output_type == 'web'):
      html_file.write ('<tr><td>scale-resolved true imbalances<br>\n')
    NumScales = (len(TRUTH_RMS) - 12) / 2
    print 'There are ', NumScales, ' different scale filters'
    lineone = []
    lineone.append('solid')
    lineone.append('dashed')
    lineone.append('dotted')

    # Plot the geostrophic imbalance
    fig, ax1 = plt.subplots()
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax2      = ax1.twinx()
    #ax1.set_ylim([0.0,0.9])
    for scale in range (NumScales):
      ax1.plot(time[:], TRUTH_RMS[12+2*scale][:], color='blue', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale] + 'km')
    ax1.legend(loc='upper left')
    #ax1.legend(loc=(0.3, legypos[fileno]))
    ax1.set_xlabel('time (hours)', color='black', fontsize=16)
    ax1.set_ylabel('Geostrophic imbalance', color='blue', fontsize=16)
    for tl in ax1.get_yticklabels():
      tl.set_color('blue')

    # Plot the hydrostatic imbalance
    #ax2.set_ylim([0.0,0.12])
    for scale in range (NumScales):
      ax2.plot(time[:], TRUTH_RMS[12+2*scale+1][:], color='red', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale+1] + 'km')
    ax2.legend(loc='lower right')
    #ax2.legend(loc=(0.6, legypos[fileno]))
    ax2.set_ylabel('Hydrostatic imbalance', color='red', fontsize=16)
    for tl in ax2.get_yticklabels():
      tl.set_color('red')
    plt.title('Geo and hydro imbalance truth', color='black', fontsize=16)

    #plt.show()
    graphics_file_name = 'Plots/ImbalScale_truth' + filesuffix
    plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
    plt.close()
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')
      html_file.write ('</tr>\n')


  # ===== Plot the data and the free model run =====
  # Plot data only if code is 'bg' or 'anal'
  if (code == 'bg') or (code == 'anal'):
    # All wind components, data
    # --------------------------
    print 'Data winds'
    if (output_type == 'web'):
      html_file.write ('<tr><td>Wind component ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m/s)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[0][:],  linewidth=1, ls='solid', color='red', label='u RMS')
    ax.plot(time[:], DATA_RMS[1][:],  linewidth=1, ls='solid', color='blue', label='v RMS')
    ax.plot(time[:], DATA_RMS[2][:],  linewidth=1, ls='solid', color='green', label='w RMS')
    ax.plot(time[:], DATA_mean[0][:], linewidth=2, ls='dotted', color='red', label='u mean')
    ax.plot(time[:], DATA_mean[1][:], linewidth=2, ls='dotted', color='blue', label='v mean')
    ax.plot(time[:], DATA_mean[2][:], linewidth=2, ls='dotted', color='green', label='w mean')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[0][:],  linewidth=1, ls='dashed', color='red', label='u RMS (free)')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[1][:],  linewidth=1, ls='dashed', color='blue', label='v RMS (free)')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[2][:],  linewidth=1, ls='dashed', color='green', label='w RMS (free)')
    ax.legend(loc='best') #(loc='upper right')
    plt.title('RMS and mean ' + code + ' values of wind components')
    graphics_file_name = 'Plots/Wind_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # rp, data
    # -------------------
    print 'Data rp'
    if (output_type == 'web'):
      html_file.write ('<tr><td>rp ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[3][:], linewidth=1, ls='solid', color='red', label='RMS')
    ax.plot(time[:], DATA_mean[3][:], linewidth=2, ls='dotted', color='red', label='mean')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[3][:],  linewidth=1, ls='dashed', color='red', label='RMS (free)')
    plt.title('RMS and mean ' + code + ' values of rp')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/rp_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # bp, data
    # -------------------
    print 'Data bp'
    if (output_type == 'web'):
      html_file.write ('<tr><td>bp ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m/s2)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[4][:], linewidth=1, ls='solid', color='blue', label='RMS')
    ax.plot(time[:], DATA_mean[4][:], linewidth=2, ls='dotted', color='blue', label='mean')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[4][:],  linewidth=1, ls='dashed', color='blue', label='RMS (free)')
    plt.title('RMS and mean ' + code + ' values of buoyancy')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/bp_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # tracer, data
    # -------------------
    print 'Data tracer'
    if (output_type == 'web'):
      html_file.write ('<tr><td>tracer ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[5][:], linewidth=1, ls='solid', color='blue', label='RMS')
    ax.plot(time[:], DATA_mean[5][:], linewidth=2, ls='dotted', color='blue', label='mean')
    ax.plot(cycle_bound_times[:-1], FREE_RMS[5][:],  linewidth=1, ls='dashed', color='blue', label='RMS (free)')
    plt.title('RMS and mean ' + code + ' values of tracer')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/tracer_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # imbalances, data
    # -------------------
    print 'Data imbalances'
    if (output_type == 'web'):
      html_file.write ('<tr><td>imbalances ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[6][:], linewidth=1, ls='solid', color='blue', label='Geo imbal RMS')
    ax.plot(time[:], DATA_RMS[7][:], linewidth=1, ls='solid', color='orange', label='Hyd imbal RMS')
    ax.plot(time[:], DATA_mean[6][:], linewidth=2, ls='dotted', color='blue', label='Hyd imbal mean')
    ax.plot(time[:], DATA_mean[7][:], linewidth=2, ls='dotted', color='orange', label='Geo imbal mean')
    plt.title('RMS and mean ' + code + ' values of imbalances')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/imbal_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # energy, data
    # -------------------
    print 'Data energy'
    if (output_type == 'web'):
      html_file.write ('<tr><td>energy ' + code + '<br>\n')
    fig, ax = plt.subplots()
    ax.set_xlabel('time (h)')
    ax.set_ylabel('RMS and mean (m2/s2)')
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax.plot(time[:], DATA_RMS[8][:], linewidth=2, ls='solid', color='red', label='KE RMS')
    ax.plot(time[:], DATA_RMS[9][:], linewidth=2, ls='solid', color='blue', label='BE RMS')
    ax.plot(time[:], DATA_RMS[10][:], linewidth=2, ls='solid', color='green', label='EE RMS')
    ax.plot(time[:], DATA_RMS[11][:], linewidth=2, ls='solid', color='black', label='TE RMS')
    ax.plot(time[:], DATA_mean[8][:], linewidth=2, ls='dotted', color='red', label='KE mean')
    ax.plot(time[:], DATA_mean[9][:], linewidth=2, ls='dotted', color='blue', label='BE mean')
    ax.plot(time[:], DATA_mean[10][:], linewidth=2, ls='dotted', color='green', label='EE mean')
    ax.plot(time[:], DATA_mean[11][:], linewidth=2, ls='dotted', color='black', label='TE mean')
    plt.title('RMS and mean ' + code + ' values of energy')
    ax.legend(loc='best') #(loc='upper right')
    graphics_file_name = 'Plots/energy_' + code + filesuffix
    #plt.show()             # Comment out to not display plot on the screen
    plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
    plt.close('all')
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

    # Scale-resolved balances, data (RMS only)
    # -------------------
    print 'Data scale resolved balances'
    if (output_type == 'web'):
      html_file.write ('<tr><td>scale-resolved ' + code + ' imbalances<br>\n')
    NumScales = (len(DATA_RMS) - 12) / 2
    print 'There are ', NumScales, ' different scale filters'
    lineone = []
    lineone.append('solid')
    lineone.append('dashed')
    lineone.append('dotted')

    # Plot the geostrophic imbalance
    fig, ax1 = plt.subplots()
    fig.set_size_inches(14.0, 7.0)
    for cycle in cycle_bound_times:
      plt.axvline(x=cycle, color='yellow')
    ax2      = ax1.twinx()
    #ax1.set_ylim([0.0,0.9])
    for scale in range (NumScales):
      ax1.plot(time[:], DATA_RMS[12+2*scale][:], color='blue', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale] + 'km')
    ax1.legend(loc='upper left')
    #ax1.legend(loc=(0.3, legypos[fileno]))
    ax1.set_xlabel('time (hours)', color='black', fontsize=16)
    ax1.set_ylabel('Geostrophic imbalance', color='blue', fontsize=16)
    for tl in ax1.get_yticklabels():
      tl.set_color('blue')

    # Plot the hydrostatic imbalance
    #ax2.set_ylim([0.0,0.12])
    for scale in range (NumScales):
      ax2.plot(time[:], DATA_RMS[12+2*scale+1][:], color='red', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale+1] + 'km')
    ax2.legend(loc='lower right')
    #ax2.legend(loc=(0.6, legypos[fileno]))
    ax2.set_ylabel('Hydrostatic imbalance', color='red', fontsize=16)
    for tl in ax2.get_yticklabels():
      tl.set_color('red')
    plt.title('Geo and hydro imbalance ' + code, color='black', fontsize=16)

    #plt.show()
    graphics_file_name = 'Plots/ImbalScale_' + code + filesuffix
    plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
    plt.close()
    if (output_type == 'web'):
      html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')
      html_file.write ('</tr>\n')



  # ===== Plot the data errors =====
  # All wind components, errors
  # --------------------------
  print 'Wind error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>Wind component ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error (m/s)')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[0][:],  linewidth=1, ls='solid', color='red', label='u RMS err')
  ax.plot(time[:], DATA_RMSE[1][:],  linewidth=1, ls='solid', color='blue', label='v RMS err')
  ax.plot(time[:], DATA_RMSE[2][:],  linewidth=1, ls='solid', color='green', label='w RMS err')
  ax.plot(time[:], DATA_err[0][:], linewidth=2, ls='dotted', color='red', label='u mean err')
  ax.plot(time[:], DATA_err[1][:], linewidth=2, ls='dotted', color='blue', label='v mean err')
  ax.plot(time[:], DATA_err[2][:], linewidth=2, ls='dotted', color='green', label='w mean err')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[0][:],  linewidth=1, ls='dashed', color='red', label='u RMS err (free)')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[1][:],  linewidth=1, ls='dashed', color='blue', label='v RMS err (free)')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[2][:],  linewidth=1, ls='dashed', color='green', label='w RMS err (free)')
  ax.legend(loc='best') #(loc='upper right')
  plt.title('RMS and mean ' + code + ' error values of wind components')
  graphics_file_name = 'Plots/Wind_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # rp, errors
  # -------------------
  print 'rp error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>rp ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[3][:], linewidth=1, ls='solid', color='red', label='RMS err')
  ax.plot(time[:], DATA_err[3][:], linewidth=2, ls='dotted', color='red', label='mean err')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[3][:],  linewidth=1, ls='dashed', color='red', label='RMS err (free)')
  plt.title('RMS and mean ' + code + ' error values of rp')
  ax.legend(loc='best') #(loc='upper right')
  graphics_file_name = 'Plots/rp_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # bp, errors
  # -------------------
  print 'bp error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>bp ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error(m/s2)')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[4][:], linewidth=1, ls='solid', color='blue', label='RMS err')
  ax.plot(time[:], DATA_err[4][:], linewidth=2, ls='dotted', color='blue', label='mean err')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[4][:],  linewidth=1, ls='dashed', color='blue', label='RMS err (free)')
  plt.title('RMS and mean ' + code + ' error values of buoyancy')
  ax.legend(loc='best') #(loc='upper right')
  graphics_file_name = 'Plots/bp_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # tracer, errors
  # -------------------
  print 'tracer error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>tracer ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[5][:], linewidth=1, ls='solid', color='blue', label='RMS err')
  ax.plot(time[:], DATA_err[5][:], linewidth=2, ls='dotted', color='blue', label='mean err')
  ax.plot(cycle_bound_times[:-1], FREE_RMSE[5][:],  linewidth=1, ls='dashed', color='blue', label='RMS err (free)')
  plt.title('RMS and mean ' + code + ' error values of tracer')
  ax.legend(loc='best') #(loc='upper right')
  graphics_file_name = 'Plots/tracer_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # imbalances, errors
  # -------------------
  print 'imbalances error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>imbalances ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[6][:], linewidth=1, ls='solid', color='blue', label='Geo imbal RMS err')
  ax.plot(time[:], DATA_RMSE[7][:], linewidth=1, ls='solid', color='orange', label='Hyd imbal RMS err')
  ax.plot(time[:], DATA_err[6][:], linewidth=2, ls='dotted', color='blue', label='Hyd imbal mean err')
  ax.plot(time[:], DATA_err[7][:], linewidth=2, ls='dotted', color='orange', label='Geo imbal mean err')
  plt.title('RMS and mean ' + code + ' error values of imbalances')
  ax.legend(loc='best') #(loc='upper right')
  graphics_file_name = 'Plots/imbal_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # energy, errors
  # -------------------
  print 'energy error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>energy ' + code + ' errors<br>\n')
  fig, ax = plt.subplots()
  ax.set_xlabel('time (h)')
  ax.set_ylabel('RMS and mean error (m2/s2)')
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax.plot(time[:], DATA_RMSE[8][:], linewidth=2, ls='solid', color='red', label='KE RMS err')
  ax.plot(time[:], DATA_RMSE[9][:], linewidth=2, ls='solid', color='blue', label='BE RMS err')
  ax.plot(time[:], DATA_RMSE[10][:], linewidth=2, ls='solid', color='green', label='EE RMS err')
  ax.plot(time[:], DATA_RMSE[11][:], linewidth=2, ls='solid', color='black', label='TE RMS err')
  ax.plot(time[:], DATA_err[8][:], linewidth=2, ls='dotted', color='red', label='KE mean err')
  ax.plot(time[:], DATA_err[9][:], linewidth=2, ls='dotted', color='blue', label='BE mean err')
  ax.plot(time[:], DATA_err[10][:], linewidth=2, ls='dotted', color='green', label='EE mean err')
  ax.plot(time[:], DATA_err[11][:], linewidth=2, ls='dotted', color='black', label='TE mean err')
  plt.title('RMS and mean ' + code + ' error values of energy')
  ax.legend(loc='best') #(loc='upper right')
  graphics_file_name = 'Plots/energy_' + code + '_err' + filesuffix
  #plt.show()             # Comment out to not display plot on the screen
  plt.savefig(plot_dir + '/' + graphics_file_name)  # Comment out to not save the file
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')

  # Scale-resolved balances, errors (RMS only)
  # -------------------
  print 'scaled resolved balance error'
  if (output_type == 'web'):
    html_file.write ('<tr><td>scale-resolved ' + code + ' imbalance errors<br>\n')
  NumScales = (len(DATA_RMSE) - 12) / 2
  print 'There are ', NumScales, ' different scale filters'
  lineone = []
  lineone.append('solid')
  lineone.append('dashed')
  lineone.append('dotted')

  # Plot the geostrophic imbalance
  fig, ax1 = plt.subplots()
  fig.set_size_inches(14.0, 7.0)
  for cycle in cycle_bound_times:
    plt.axvline(x=cycle, color='yellow')
  ax2      = ax1.twinx()
  #ax1.set_ylim([0.0,0.9])
  for scale in range (NumScales):
    ax1.plot(time[:], DATA_RMSE[12+2*scale][:], color='blue', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale] + 'km')
  ax1.legend(loc='upper left')
  #ax1.legend(loc=(0.3, legypos[fileno]))
  ax1.set_xlabel('time (hours)', color='black', fontsize=16)
  ax1.set_ylabel('Geostrophic imbalance error', color='blue', fontsize=16)
  for tl in ax1.get_yticklabels():
    tl.set_color('blue')
    
  # Plot the hydrostatic imbalance
  #ax2.set_ylim([0.0,0.12])
  for scale in range (NumScales):
    ax2.plot(time[:], DATA_RMSE[12+2*scale+1][:], color='red', linewidth='2', ls=lineone[scale], label='>' + quantities[12+2*scale+1] + 'km')
  ax2.legend(loc='lower right')
  #ax2.legend(loc=(0.6, legypos[fileno]))
  ax2.set_ylabel('Hydrostatic imbalance error', color='red', fontsize=16)
  for tl in ax2.get_yticklabels():
    tl.set_color('red')
  plt.title('Geo and hydro imbalance ' + code + 'error', color='black', fontsize=16)

  #plt.show()
  graphics_file_name = 'Plots/ImbalScale_' + code + '_err' + filesuffix
  plt.savefig(plot_dir + '/' + graphics_file_name, bbox_inches='tight')
  plt.close()
  if (output_type == 'web'):
    html_file.write ('<img src=' + graphics_file_name + ' width=800></td>\n')
    html_file.write ('</tr>\n')#
    html_file.write ('</table>')

  return





# ===================================================================
# ===================================================================
def make_histogram_ob (nbins, obs, back_obs, anal_obs, true_obs, quantity, output_type, plot_dir, html_file):
  # Subroutine to generate histogram data of observations - background observations,
  #                                          observations - analysis observations,
  #                                          observations - truth,
  #                                          background observations - truth,
  #                                          analysis observations - truth

  if (output_type == 'web'):
    filesuffix = '.png'
    html_file.write ('<h4>Observation differences for ' + quantity + '</h4>\n')
    html_file.write ('<table>\n')
    html_file.write ('<tr>\n')
  else:
    filesuffix = '.eps'

  # Number of observations
  Nobs = len(obs)

  # Convert the inputs to numpy arrays
  obs      = np.asarray(obs)
  back_obs = np.asarray(back_obs)
  anal_obs = np.asarray(anal_obs)
  true_obs = np.asarray(true_obs)

  # Compute the differences
  omb      = obs      - back_obs
  oma      = obs      - anal_obs
  omt      = obs      - true_obs
  bmt      = back_obs - true_obs
  amt      = anal_obs - true_obs

  # Compute the means
  omb_mean = np.mean(omb)
  oma_mean = np.mean(oma)
  omt_mean = np.mean(omt)
  bmt_mean = np.mean(bmt)
  amt_mean = np.mean(amt)

  # Compute the standard deviations
  omb_std = np.std(omb)
  oma_std = np.std(oma)
  omt_std = np.std(omt)
  bmt_std = np.std(bmt)
  amt_std = np.std(amt)

  #print '----- Obs - bg'
  #print omb
  #print '----- Obs - an'
  #print oma
  #print '----- Obs - tr'
  #print omt
  #print '----- Bg - tr'
  #print bmt
  #print '----- An - tr'
  #print amt

  # Find the maximum absolute value of these differences
  maxval = 0.0
  maxval = max([maxval, max(abs(omb))])
  maxval = max([maxval, max(abs(oma))])
  maxval = max([maxval, max(abs(omt))])
  maxval = max([maxval, max(abs(bmt))])
  maxval = max([maxval, max(abs(amt))])
  maxval = flexi_round(maxval)

  # Generate the bins
  bins   = np.linspace(-1.0*maxval, maxval, nbins+1)
  print 'The following bins have been generated'
  print bins

  # Generate the histograms (not needed as the work is done inside ax.hist below)
  #print '----- Generating histogram for omb'
  #hist_omb = gen_histogram (bins, omb)
  #print hist_omb
  #print '----- Generating histogram for oma'
  #hist_oma = gen_histogram (bins, oma)
  #print hist_oma
  #print '----- Generating histogram for omt'
  #hist_omt = gen_histogram (bins, omt)
  #print hist_omt
  #print '----- Generating histogram for bmt'
  #hist_bmt = gen_histogram (bins, bmt)
  #print hist_bmt
  #print '----- Generating histogram for amt'
  #hist_amt = gen_histogram (bins, amt)
  #print hist_amt

  # Plot histograms

  # Observations - background
  fig, ax = plt.subplots()
  ax.hist(omb, bins=bins, histtype='bar', facecolor='b')
  ax.set_title('O-B for ' + quantity + '\nmean=' + str.format('%.5f' % omb_mean) + ', stddev=' + str.format('%.5f' % omb_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-b', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = 'Plots/' + quantity + '_o-b' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td>Observations - background<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')


  # Observations - analysis
  fig, ax = plt.subplots()
  ax.hist(oma, bins=bins, histtype='bar', facecolor='r')
  ax.set_title('O-A for ' + quantity + '\nmean=' + str.format('%.5f' % oma_mean) + ', stddev=' + str.format('%.5f' % oma_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-a', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = 'Plots/' + quantity + '_o-a' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td>Observations - analysis<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Observations - truth
  fig, ax = plt.subplots()
  ax.hist(omt, bins=bins, histtype='bar', facecolor='g')
  ax.set_title('O-T for ' + quantity + '\nmean=' + str.format('%.5f' % omt_mean) + ', stddev=' + str.format('%.5f' % omt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('o-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = 'Plots/' + quantity + '_o-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td>Observations - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Background - truth
  fig, ax = plt.subplots()
  ax.hist(bmt, bins=bins, histtype='bar', facecolor='cyan')
  ax.set_title('B-T for ' + quantity + '\nmean=' + str.format('%.5f' % bmt_mean) + ', stddev=' + str.format('%.5f' % bmt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('b-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = 'Plots/' + quantity + '_b-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td>Background - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')

  # Analysis - truth
  fig, ax = plt.subplots()
  ax.hist(amt, bins=bins, histtype='bar', facecolor='purple')
  ax.set_title('A-T for ' + quantity + '\nmean=' + str.format('%.5f' % amt_mean) + ', stddev=' + str.format('%.5f' % amt_std) + ', Nobs=' + str(Nobs), fontsize=16)
  ax.set_xlabel('a-t', fontsize=16)
  ax.set_ylabel('frequency', fontsize=16)
  figfilename = 'Plots/' + quantity + '_a-t' + filesuffix
  plt.savefig(plot_dir + '/' + figfilename, bbox_inches='tight')
  plt.close('all')
  if (output_type == 'web'):
    html_file.write ('<td>Analysis - truth<br>\n')
    html_file.write ('<img src=' + figfilename  + ' width=300></td>\n')
    html_file.write ('</tr>\n')
    html_file.write ('</table>\n')
  return


# ===================================================================
# ===================================================================
def find_index (bounds, value):
  # Find the lower index of a value

  arraylen = len(bounds) - 1
  here     = 0
  vhere    = bounds[here]
  while (vhere <= value) and (here < arraylen):
    here += 1
    if here < arraylen:
      vhere = bounds[here]

  if here <= arraylen:
    here -= 1

#  if here > arraylen:
#    print value, ' is beyond ', bounds[arraylen]
#  else:
#    here -= 1
#    print value, ' is between ', bounds[here], ' and ', bounds[here+1]

  return here


# ===================================================================
# ===================================================================
def gen_histogram (bins, values):
  # To generate a histogram

  length = len(bins)
  hist   = np.repeat(0, length)

  for value in values:
    place = find_index(bins, value)
    if place < length:
      hist[place] += 1
    else:
      print 'Out of bounds'
  return hist

# ===================================================================
# ===================================================================
