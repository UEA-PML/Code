#!/usr/bin/env python
import sys
import os
import datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib import rc
from datetime import datetime
from datetime import date, timedelta
from numpy import copy, nanmean, array
import pylab 
import math
import argparse
import pandas as pd
from scipy import signal
from read_instruments import *
from EC_functions import *
import argparse

rc('text',usetex=True)

def main():

   parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter
   )

   parser.add_argument("-c", "--CR6",  action='store_true',  default=False, help="Use CR6 logged LPMS motion data rather than SysTron Donner.")
   parser.add_argument("-l", "--LPMS",  action='store_true',  default=False, help="Use LPMS motion data rather than SysTron Donner.")
   parser.add_argument("-o", "--ofile", type=str, default='', help="Output filename")
   parser.add_argument("-nl", "--no_licor", action='store_true',  default=False, help="Use if no LICOR data available")
   parser.add_argument("-u", "--underway", action='store_true', default=False, help="Use ship's underway data if available")
   parser.add_argument("-L0", "--L0", type=str, default='', help="L0 Output filestem")
   parser.add_argument("-g", "--gill", action='store_true',  default=False, help="Use Gill anemometer data.")
   parser.add_argument("-aco", "--ACO", action='store_true',  default=False, help="Navigational data in ACO format.")
   args = parser.parse_args()

   ######## Definition of fixed parameters ##########
   deg2rad = math.pi/180.    # conversion of degrees to radians
   rad2deg = 180./math.pi    # conversion of radians to degrees
   g  = 9.80665              # gravitational acceleration (m/s2) - conventional standard value
   knots2ms = 0.51444444     # conversion of knots to m/s
   deg2K = 273.15            # conversion of temperature degC to Kelvin (K = degC + 273.15)
   Cp_air = 1004.            # heat capacity dry air (J/kg/K)
   Rv = 287.058              # specific gas constant for dry air (J/kg/K)
   P_std = 101.325e+3        # standard atmospheric pressure (Pa)
   Ca = 2.16679e-3           # Absolute humidity calculation (kg.K/J)
   g_kg = 1e+3               # grams (g) in a kilogram (kg)
   molCO2 = 44.0095          # molecular weight of CO2
   molCH4 = 16.04246         # molecular weight of CH4
   ppm = 1e+6                # parts per million 
   mmol = 1e+3               # milli mols to mol
   umol = 1e+6               # micro mols to mol
   nmol = 1e+9               # nano mols to mol
   Licor_lag_s = 0.2         # Lag in Licor data due to processing (0.2s)
   Picarro_lag_s = 3.3       # Hardwired value to account for transit time between inlet and instrument and 
                             # calculated from data 12 Oct - 28 Oct 2018
                             # This is replaced in the code when there is a puff test
   RHO_AIR = 1.2             # Density of air in kg/m3
   R_star = 8.3143           # Universal gas constant (R*) (J/K/mol)
   T_std = 293.15            # standard temperature (K)
   T_std_C = T_std - 273.15  # standard temperature (deg C)
   M_w = 18.016              # molecular weight of water
   seconds_per_day = 24.*60.*60. # seconds in 1 day (s/day)

   # Calculate the latent heat of vaporisation (as this has a temperature dependency)
   Lvap = (2.501 - 0.00237*T_std_C)*1e+3   # J/g
   # ideal gas at standard pressure and temperature
   PRT = P_std/(R_star*T_std)

   # Define position vector for distance between MotionPak and sonic volume.
   # Units are metres, x positive forward, y positive to port, z positive up.
   # These values are for the system on the JCR
   Rx = 0.17   #Metek Sonic is 17 cm in front of MotionPak
   Ry = 0.  
   Rz = 0.725  #Metek Sonic is 72.5 cm above MotionPak
   alpha = -0.3 #Metek Sonic rotated 0.3 degrees towards port side relative to the motion sensor
   beta = 1.2  #Metek Sonic rotated 1.2 degrees towards the stern relative to the motion sensor

   fsam = 10   # All data to be interpolated onto a 10 Hz grid
   if (fsam == 10):
      freq_str = '100L' # interpreted as 100 ms (0.1 s)

   ##########
   
   #### PART 1 #####
   # Reading in the various sensor data streams
   #   i. Wind data (function - read_wind() for Sonic; read_gill() for Gill)
   #  ii. LPMS data (function - read_LPMS())
   # iii. Sysdon data (function - read_sysdon())
   #  iv. Picarro data (function - read_picarro())
   #   v. Licor data (function - read_licor())
   #  vi. seatex data (functions - read_seatex_gll() & read_seatex_vtg())
   # vii. gyro data (function - read_gyro())
   # 
   # Following data input all data to be used in onward processing is put onto a
   # consistent grid (resampled at rate 'fsam' defined above - typically 10Hz)
   
   # i. Wind data
   print "Reading in wind data"
   if args.gill:
      print " Using Gill anemometer (1 Hz)"
      # approximate positions of the Gill
      Rx = 0.10
      Ry = -0.6
      Rz = 0.8
      alpha = -30.0
      beta = 0.0
      u_raw, v_raw, w_raw, t_degC, datetimeList_wind = read_gill()
      polarity = 0 # no need to reverse polarity for the Gill
      # correct Gill data for angular offset of installation
      u_ms, v_ms, w_ms = angular_offset(u_raw, v_raw, w_raw, alpha*deg2rad, beta*deg2rad, polarity)
   else:
      print " Using Metek Sonic anemometer (10 Hz)"
      u_raw, v_raw, w_raw, t_degC, datetimeList_wind = read_wind()
      polarity = 1 # need to reverse polarity for Metek
      # correct Metek data for angular offset of installation
      u_ms, v_ms, w_ms = angular_offset(u_raw, v_raw, w_raw, alpha*deg2rad, beta*deg2rad, polarity)

   wind_start_time = datetimeList_wind[0]
   wind_end_time = datetimeList_wind[len(datetimeList_wind)-1]
   print "Time range: ", wind_start_time,"->", wind_end_time
   
   # Interpolate all wind variables onto a common time grid defined by fsam
   wind_u_ms_ts = pd.Series(u_ms, index = datetimeList_wind)
   wind_u_ms_ts = wind_u_ms_ts.resample(freq_str).mean().ffill()
   wind_v_ms_ts = pd.Series(v_ms, index = datetimeList_wind)
   wind_v_ms_ts = wind_v_ms_ts.resample(freq_str).mean().ffill()
   wind_w_ms_ts = pd.Series(w_ms, index = datetimeList_wind)
   wind_w_ms_ts = wind_w_ms_ts.resample(freq_str).mean().ffill()
   t_degC_ts = pd.Series(t_degC, index = datetimeList_wind)
   t_degC_ts = t_degC_ts.resample(freq_str).mean().ffill()
   
   # ii. LPMS data
   if args.LPMS:
      print "Reading in LPMS data"
      lrotx_degs, lroty_degs, lrotz_degs, laccelx_g, laccely_g, laccelz_g, datetimeList_LPMS = read_LPMS()
      LPMS_start_time = datetimeList_LPMS[0]
      LPMS_end_time = datetimeList_LPMS[len(datetimeList_LPMS)-1]
      print "Time range: ", LPMS_start_time,"->", LPMS_end_time
   
      # Interpolate all LPMS variables onto a common time grid defined by fsam
      rotx_degs_ts = pd.Series(lrotx_degs, index = datetimeList_LPMS)
      rotx_degs_ts = rotx_degs_ts.resample(freq_str).mean().ffill()
      roty_degs_ts = pd.Series(lroty_degs, index = datetimeList_LPMS)
      roty_degs_ts = roty_degs_ts.resample(freq_str).mean().ffill()
      rotz_degs_ts = pd.Series(lrotz_degs, index = datetimeList_LPMS)
      rotz_degs_ts = rotz_degs_ts.resample(freq_str).mean().ffill()
      accelx_g_ts = pd.Series(laccelx_g, index = datetimeList_LPMS)
      accelx_g_ts = accelx_g_ts.resample(freq_str).mean().ffill()
      accely_g_ts = pd.Series(laccely_g, index = datetimeList_LPMS)
      accely_g_ts = accely_g_ts.resample(freq_str).mean().ffill()
      accelz_g_ts = pd.Series(laccelz_g, index = datetimeList_LPMS)
      accelz_g_ts = accelz_g_ts.resample(freq_str).mean().ffill()

   if args.CR6:
      print "Reading in CR6 data"
      u_raw, v_raw, w_raw, t_degC, lrotx_degs, lroty_degs, lrotz_degs, laccelx_g, laccely_g, laccelz_g, fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList_CR6 = read_CR6()

      CR6_start_time = datetimeList_CR6[0]
      CR6_end_time = datetimeList_CR6[len(datetimeList_CR6)-1]
      print "Time range: ", CR6_start_time,"->", CR6_end_time

      polarity = 1 # need to reverse polarity for Metek
      # correct Metek data for angular offset of installation
      u_ms, v_ms, w_ms = angular_offset(u_raw, v_raw, w_raw, alpha*deg2rad, beta*deg2rad, polarity)

      # Interpolate all LPMS variables onto a common time grid defined by fsam
      rotx_degs_ts = pd.Series(lrotx_degs, index = datetimeList_CR6)
      rotx_degs_ts = rotx_degs_ts.resample(freq_str).mean().ffill()
      roty_degs_ts = pd.Series(lroty_degs, index = datetimeList_CR6)
      roty_degs_ts = roty_degs_ts.resample(freq_str).mean().ffill()
      rotz_degs_ts = pd.Series(lrotz_degs, index = datetimeList_CR6)
      rotz_degs_ts = rotz_degs_ts.resample(freq_str).mean().ffill()
      accelx_g_ts = pd.Series(laccelx_g, index = datetimeList_CR6)
      accelx_g_ts = accelx_g_ts.resample(freq_str).mean().ffill()
      accely_g_ts = pd.Series(laccely_g, index = datetimeList_CR6)
      accely_g_ts = accely_g_ts.resample(freq_str).mean().ffill()
      accelz_g_ts = pd.Series(laccelz_g, index = datetimeList_CR6)
      accelz_g_ts = accelz_g_ts.resample(freq_str).mean().ffill()
      
      # Interpolate all wind variables onto a common time grid defined by fsam
      CR6_wind_u_ms_ts = pd.Series(u_ms, index = datetimeList_CR6)
      CR6_wind_u_ms_ts = CR6_wind_u_ms_ts.resample(freq_str).mean().ffill()
      CR6_wind_v_ms_ts = pd.Series(v_ms, index = datetimeList_CR6)
      CR6_wind_v_ms_ts = CR6_wind_v_ms_ts.resample(freq_str).mean().ffill()
      CR6_wind_w_ms_ts = pd.Series(w_ms, index = datetimeList_CR6)
      CR6_wind_w_ms_ts = CR6_wind_w_ms_ts.resample(freq_str).mean().ffill()
      CR6_t_degC_ts = pd.Series(t_degC, index = datetimeList_CR6)
      CR6_t_degC_ts = CR6_t_degC_ts.resample(freq_str).mean().ffill()

   # iii. SysTron Donner data
   if not(args.LPMS or args.CR6):
      print "Reading in SysDon Data"
      srotx_degs, sroty_degs, srotz_degs, saccelx_g, saccely_g, saccelz_g, datetimeList_sysdon = read_sysdon()
      sysdon_start_time = datetimeList_sysdon[0]
      sysdon_end_time = datetimeList_sysdon[len(datetimeList_sysdon)-1]
      print "Time range: ", sysdon_start_time,"->", sysdon_end_time
      
      # Interpolate all sysdon variables onto a common time grid defined by fsam
      rotx_degs_ts = pd.Series(srotx_degs, index = datetimeList_sysdon)
      rotx_degs_ts = rotx_degs_ts.resample(freq_str).mean().ffill()
      roty_degs_ts = pd.Series(sroty_degs, index = datetimeList_sysdon)
      roty_degs_ts = roty_degs_ts.resample(freq_str).mean().ffill()
      rotz_degs_ts = pd.Series(srotz_degs, index = datetimeList_sysdon)
      rotz_degs_ts = rotz_degs_ts.resample(freq_str).mean().ffill()
      accelx_g_ts = pd.Series(saccelx_g, index = datetimeList_sysdon)
      accelx_g_ts = accelx_g_ts.resample(freq_str).mean().ffill()
      accely_g_ts = pd.Series(saccely_g, index = datetimeList_sysdon)
      accely_g_ts = accely_g_ts.resample(freq_str).mean().ffill()
      accelz_g_ts = pd.Series(saccelz_g, index = datetimeList_sysdon)
      accelz_g_ts = accelz_g_ts.resample(freq_str).mean().ffill()
   
   # iv. Licor data
   if not(args.no_licor):
      print "Reading in Licor data"
      fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList_licor = read_licor()
      # add the time offset (for processing data)
      timedelta_Licor = timedelta(seconds=Licor_lag_s)
      for row in range(len(datetimeList_licor)):
         datetimeList_licor[row] = datetimeList_licor[row] - timedelta_Licor
      licor_start_time = datetimeList_licor[0]
      licor_end_time = datetimeList_licor[len(datetimeList_licor)-1]
      print "Time range: ", licor_start_time,"->",licor_end_time 
   else:
      print "No Licor data being used - dummy arrays created"
      fco2_density = np.zeros(len(u_ms))*0.
      fh2o_density = np.zeros(len(u_ms))*0.
      fco2_ppm = np.zeros(len(u_ms))*0.
      fh2o_ppm = np.zeros(len(u_ms))*0.
      ftemp = np.zeros(len(u_ms))*0.
      fpressure = np.zeros(len(u_ms))*0.
      fdiag = np.zeros(len(u_ms))*0.
      datetimeList_licor = datetimeList_wind
   
   # Interpolate all licor variables onto a common time grid defined by fsam
   licor_fco2_density_ts = pd.Series(fco2_density, index = datetimeList_licor)
   licor_fco2_density_ts = licor_fco2_density_ts.resample(freq_str).mean().ffill()
   licor_fh2o_density_ts = pd.Series(fh2o_density, index = datetimeList_licor)
   licor_fh2o_density_ts = licor_fh2o_density_ts.resample(freq_str).mean().ffill()
   licor_fco2_ppm_ts = pd.Series(fco2_ppm, index = datetimeList_licor)
   licor_fco2_ppm_ts = licor_fco2_ppm_ts.resample(freq_str).mean().ffill()
   licor_fh2o_ppm_ts = pd.Series(fh2o_ppm, index = datetimeList_licor)
   licor_fh2o_ppm_ts = licor_fh2o_ppm_ts.resample(freq_str).mean().ffill()
   licor_ftemp_ts = pd.Series(ftemp, index = datetimeList_licor)
   licor_ftemp_ts = licor_ftemp_ts.resample(freq_str).mean().ffill()
   licor_fpressure_ts = pd.Series(fpressure, index = datetimeList_licor)
   licor_fpressure_ts = licor_fpressure_ts.resample(freq_str).mean().ffill()
   licor_fdiag_ts = pd.Series(fdiag, index = datetimeList_licor)
   licor_fdiag_ts = licor_fdiag_ts.resample(freq_str).mean().ffill()

   # v. seatex data
   # Reading in the ancillary data for position, heading and speed
   if (args.ACO):
      print "Reading in seatex (gga) data for position"
      lat, lon, datetimeList_seatex_gga = read_seatex_gga()
      gga_start_time = datetimeList_seatex_gga[0]
      gga_end_time = datetimeList_seatex_gga[len(datetimeList_seatex_gga)-1]
      print "Time range: ", gga_start_time,"->", gga_end_time
   else:   
      print "Reading in seatex (gll) data for position"
      lat, lon, datetimeList_seatex_gll = read_seatex_gll()
      seatex_gll_start_time = datetimeList_seatex_gll[0]
      seatex_gll_end_time = datetimeList_seatex_gll[len(datetimeList_seatex_gll)-1]
      print "Time range: ", seatex_gll_start_time, "->", seatex_gll_end_time
   
   print "Reading in seatex (vtg) data for speed"
   # Use the gyro data rather than seatex for heading as is more accurate
   heading, SOG, datetimeList_seatex_vtg = read_seatex_vtg(args.ACO)
   seatex_vtg_start_time = datetimeList_seatex_vtg[0]
   seatex_vtg_end_time = datetimeList_seatex_vtg[len(datetimeList_seatex_vtg)-1]
   print "Time range: ", seatex_vtg_start_time, "->", seatex_vtg_end_time
   
   print "Reading in gyro data for heading"
   heading, datetimeList_gyro = read_gyro(args.ACO)
   gyro_start_time = datetimeList_gyro[0]
   gyro_end_time = datetimeList_gyro[len(datetimeList_gyro)-1]
   print "Time range: ", gyro_start_time, "->", gyro_end_time
  
   if (args.underway):
      print "Reading in underway data for met / surface parameters"
      tair, rhum, TIR, atm_pr, sst, sal, datetimeList_underway = read_underway(args.ACO)
      underway_start_time = datetimeList_underway[0]
      underway_end_time = datetimeList_underway[len(datetimeList_underway)-1]
      print "Time range: ", underway_start_time, "->", underway_end_time

      # Need to make sure that underway time-series is between the periods sampled
      # by the Picarro, rather than being an entire day's worth of data
      tair_ts = pd.Series(tair, index=datetimeList_underway)
      tair_ts = tair_ts.resample(freq_str).mean().ffill()
      tair_ts = tair_ts.loc[wind_start_time:wind_end_time]

      rhum_ts = pd.Series(rhum, index=datetimeList_underway)
      rhum_ts = rhum_ts.resample(freq_str).mean().ffill()
      rhum_ts = rhum_ts.loc[wind_start_time:wind_end_time]

      TIR_ts = pd.Series(TIR, index=datetimeList_underway)
      TIR_ts = TIR_ts.resample(freq_str).mean().ffill()
      TIR_ts = TIR_ts.loc[wind_start_time:wind_end_time]

      atm_pr_ts = pd.Series(atm_pr, index=datetimeList_underway)
      atm_pr_ts = atm_pr_ts.resample(freq_str).mean().ffill()
      atm_pr_ts = atm_pr_ts.loc[wind_start_time:wind_end_time]

      sst_ts = pd.Series(sst, index=datetimeList_underway)
      sst_ts = sst_ts.resample(freq_str).mean().ffill()
      sst_ts = sst_ts.loc[wind_start_time:wind_end_time]

      sal_ts = pd.Series(sal, index=datetimeList_underway)
      sal_ts = sal_ts.resample(freq_str).mean().ffill()
      sal_ts = sal_ts.loc[wind_start_time:wind_end_time]

   else:
      tair = np.nan
      rhum = np.nan
      TIR = np.nan
      atm_pr = np.nan
      sst = np.nan
      sal = np.nan
   
   # Interpolate all navigation variables onto a common time grid defined by fsam
   # Additionally need to make sure that time-series are between the periods sampled
   # by the Picarro, rather than being an entire day's worth of data
   if args.ACO:
      datetimeList_seatex_gll = datetimeList_seatex_gga
      
   lat_ts = pd.Series(lat, index = datetimeList_seatex_gll)
   lat_ts = lat_ts.resample(freq_str).mean().ffill()
   lat_ts = lat_ts.loc[wind_start_time:wind_end_time]

   lon_ts = pd.Series(lon, index = datetimeList_seatex_gll)
   lon_ts = lon_ts.resample(freq_str).mean().ffill()
   lon_ts = lon_ts.loc[wind_start_time:wind_end_time]

   SOG_ts = pd.Series(SOG, index = datetimeList_seatex_vtg)
   SOG_ts = SOG_ts.resample(freq_str).mean().ffill()
   SOG_ts = SOG_ts.loc[wind_start_time:wind_end_time]

   heading_ts = pd.Series(heading, index = datetimeList_gyro)
   heading_ts = heading_ts.resample(freq_str).mean().ffill()
   heading_ts = heading_ts.loc[wind_start_time:wind_end_time]

   # vi. Picarro data
   print "Reading in Picarro data"
   ch4, fco2_dry, fswitch_state, fcavity_pressure_torr, datetimeList_picarro = read_picarro()
   picarro_start_time = datetimeList_picarro[0]
   picarro_end_time = datetimeList_picarro[len(datetimeList_picarro)-1]
   print "Time range: ", picarro_start_time,"->",picarro_end_time
   
   # Interpolate all Picarro variables onto a common time grid defined by fsam
   picarro_ch4_ts = pd.Series(ch4, index = datetimeList_picarro)
   picarro_ch4_ts = picarro_ch4_ts.resample(freq_str).mean().ffill()
   picarro_fco2_dry_ts = pd.Series(fco2_dry, index = datetimeList_picarro)
   picarro_fco2_dry_ts = picarro_fco2_dry_ts.resample(freq_str).mean().ffill()
   picarro_fswitch_state = pd.Series(fswitch_state, index = datetimeList_picarro)
   picarro_fswitch_state = picarro_fswitch_state.resample(freq_str).mean().ffill()
   
   # Write out the data matrix as a binary (pickle) file
   # Need to do this here before the Picarro data are time-shifted
   L0_freq_data_matrix = pd.concat([wind_u_ms_ts,wind_v_ms_ts,wind_w_ms_ts,t_degC_ts, 
                         rotx_degs_ts, roty_degs_ts, rotz_degs_ts,
                         accelx_g_ts, accely_g_ts, accelz_g_ts,
                         picarro_ch4_ts, picarro_fco2_dry_ts, picarro_fswitch_state,
                         licor_fco2_density_ts, licor_fh2o_density_ts, licor_ftemp_ts,
                         licor_fpressure_ts, licor_fdiag_ts, 
                         lat_ts, lon_ts, heading_ts, SOG_ts], 
                         axis=1).ffill().dropna()

   # Write out this file as the L0 file
   if args.L0:
      print "Writing output L0 binary file"
      L0_filename = args.L0
      output_file = open(L0_filename,'wb')
      L0_freq_data_matrix.to_pickle(output_file)
      output_file.close()

   # find where picarro switch state is on (1)
   switch_on_index = np.where(picarro_fswitch_state == 1)
   switch_on_index = switch_on_index[0]
   
   # then find when CH4 drops below 0.004 (3 sigma) of ambient
   Picarro_response_s = np.nan
   if len(switch_on_index) > 0:
      # Calculate the mean CH4 for (5 second) period before triggering
      picarro_ch4_ts_start = np.mean(picarro_ch4_ts[switch_on_index[0]-49:switch_on_index[0]])
      # Calculate the mean CH4 for (5 second period towards end of valve fire
      picarro_ch4_ts_end = np.mean(picarro_ch4_ts[switch_on_index[len(switch_on_index)-51]:switch_on_index[len(switch_on_index)-1]])
      # Threshold for instrument response time 
      picarro_ch4_response_threshold = picarro_ch4_ts_end + (picarro_ch4_ts_start - picarro_ch4_ts_end)*(1./math.e)
      for ss in range(len(switch_on_index)):
         if (picarro_ch4_ts_start - picarro_ch4_ts[switch_on_index[ss]] > 0.004):
            print "PICARRO PUFF TEST TIME LAG(S):"
            Picarro_lag_s = float(ss)/10.
            print Picarro_lag_s
            break
      for ss in range(len(switch_on_index)):
         if (picarro_ch4_ts[switch_on_index[ss]] < picarro_ch4_response_threshold):
            print "Picarro response time(s):"
            Picarro_response_s = float(ss)/10. - Picarro_lag_s
            print Picarro_response_s
            break

   # Picarro data lags behind the ship motion
   # Either use the predetermined value or the value determined from the data when N puff activated
   timedelta_Picarro = timedelta(seconds=Picarro_lag_s)
   print "Original first data point in picarro"
   print datetimeList_picarro[0]
   for row in range(len(datetimeList_picarro)):
      datetimeList_picarro[row] = datetimeList_picarro[row] - timedelta_Picarro

   print "Shifted first data point in picarro"
   print datetimeList_picarro[0]
   # Interpolate all Picarro variables onto a common time grid defined by fsam
   picarro_ch4_ts = pd.Series(ch4, index = datetimeList_picarro)
   picarro_ch4_ts = picarro_ch4_ts.resample(freq_str).mean().ffill()
   picarro_fco2_dry_ts = pd.Series(fco2_dry, index = datetimeList_picarro)
   picarro_fco2_dry_ts = picarro_fco2_dry_ts.resample(freq_str).mean().ffill()
   picarro_fswitch_state = pd.Series(fswitch_state, index = datetimeList_picarro)
   picarro_fswitch_state = picarro_fswitch_state.resample(freq_str).mean().ffill()

   #### PART 2 ####
   # Getting data from different sources (sonic, motion, picarro, licor, navigation) onto a consistent grid
   #
   # This grid is a python *data frame*, so need to take care when extracting
   # elements from it.
   freq_data_matrix = pd.concat([wind_u_ms_ts,wind_v_ms_ts,wind_w_ms_ts,t_degC_ts, 
                      rotx_degs_ts, roty_degs_ts, rotz_degs_ts,
                      accelx_g_ts, accely_g_ts, accelz_g_ts,
                      picarro_ch4_ts, picarro_fco2_dry_ts, picarro_fswitch_state,
                      licor_fco2_density_ts, licor_fh2o_density_ts, licor_ftemp_ts,
                      licor_fpressure_ts, licor_fdiag_ts, 
                      lat_ts, lon_ts, heading_ts, SOG_ts], 
                      axis=1).ffill().dropna()

   # WindMaster Orientation already consistent with NOAA sign convention
   # U: forward is positive
   # V: to port is positive
   # W: up is positive
   U = copy(freq_data_matrix.values[:,0])
   V = copy(freq_data_matrix.values[:,1])
   W = copy(freq_data_matrix.values[:,2])
   T_K = copy(freq_data_matrix.values[:,3]) + deg2K
         
   X_rot = copy(freq_data_matrix.values[:,4])
   Y_rot = copy(freq_data_matrix.values[:,5])
   Z_rot = copy(freq_data_matrix.values[:,6])
   X_accel = copy(freq_data_matrix.values[:,7])
   Y_accel = copy(freq_data_matrix.values[:,8])
   Z_accel = copy(freq_data_matrix.values[:,9])

   # Change signs on raw motion to match NOAA sign conventions.
   # Assume x,y,z axis definitions from the NOAA script:
   #   accx  : forward is positive 
   #   accy  : to port is positive
   #   accz  : up is positive
   #   ratex : port up is positive phi
   #   ratey : bow down is positive theta
   #   ratez : bow to port is positive psi (right-handed)
   # Apply unit conversion factors to motion and heading data:
   # Angular rate must be radians/s and acceleration must be m/s^2
   ratex = -Y_rot * deg2rad
   ratey =  X_rot * deg2rad
   ratez =  Z_rot * deg2rad
   accx  = -Y_accel * g 
   accy  =  X_accel * g
   accz  =  Z_accel * g

   # Convert Speed Of Ship from knots to m/s
   smg_knots = copy(freq_data_matrix.values[:,21])
   SOS = smg_knots * knots2ms   

   # Get ship's heading in radians for use in yaw angle and true wind calculation.
   heading_deg = copy(freq_data_matrix.values[:,20])
   heading = np.unwrap(heading_deg * deg2rad) 
   
   # Extract the various gases out for the EC calculations
   CH4_picarro = copy(freq_data_matrix.values[:,10])
   CO2_picarro = copy(freq_data_matrix.values[:,11]) 
   switch_state_picarro = copy(freq_data_matrix.values[:,12])
   CO2_licor = copy(freq_data_matrix.values[:,13]) 
   H2O_licor = copy(freq_data_matrix.values[:,14])
   P_licor = copy(freq_data_matrix.values[:,16])
   sig_licor = copy(freq_data_matrix.values[:,17])
   
   # Get ship's position
   ship_lat = copy(freq_data_matrix.values[:,18])
   ship_lon = copy(freq_data_matrix.values[:,19])
   #### END OF PART 2 ####

   #### PART 3 ####
   # Blomquist motcorr function
   #
   # Readying the data matrices
   sens_disp = np.zeros(3)*0.
   sens_disp[0] = Rx
   sens_disp[1] = Ry
   sens_disp[2] = Rz
   
   # Getting the data out of the dataframe column (.axes)
   timestamps = freq_data_matrix.axes
   timestamps = timestamps[0]
   time_array_hours = timestamps.hour
   time_array_mins = timestamps.minute
   time_array_secs = timestamps.second
   time_array_msec = timestamps.microsecond

   decimal_time = np.zeros(len(time_array_hours))*0.

   #son: Nx4 array of time and sonic wind velocities.
   #    N is typically 36000 for 10 Hz data.  
   son = np.zeros((len(U),4))*0.
   for j in range(len(time_array_hours)):
      decimal_time[j] = float(time_array_hours[j]) + float(time_array_mins[j])/60. + float(time_array_secs[j])/3600. + float(time_array_msec[j])/(3600.*1e+6)
      son[j,0] = decimal_time[j]
      #son[j,1] = -1.*U[j]
      # reversal of U is presumably carried out in Blomquist corrections
      son[j,1] = U[j]
      son[j,2] = V[j]
      son[j,3] = W[j]

   #mot: Nx7 array of time, linear acceleration and angular rate variables.
   #assuming for now that the order is x,y,z
   mot = np.zeros((len(accx),7))*0.
   for j in range(len(time_array_hours)):
      mot[j,0] = decimal_time[j]
      mot[j,1] = accx[j]
      mot[j,2] = accy[j]
      mot[j,3] = accz[j]
      mot[j,4] = ratex[j]
      mot[j,5] = ratey[j]
      mot[j,6] = ratez[j]
   
   print "Performing motion correction"
   uvw_mc, acc, uvwplat, xyzplat = motcorr(son,mot,heading,sens_disp)
   uplat_Var = np.var(uvwplat[:,0])
   vplat_Var = np.var(uvwplat[:,1])
   wplat_Var = np.var(uvwplat[:,2])
   
   Zplat = xyzplat[:,2]
   Wplat = uvwplat[:,2]
   Aplat = mot[:,3]
   #### END OF PART 3 ####
   
   #### PART 4 ####
   # Calculation of wind and motion data decorrelation
   print "Performing wind and motion data decorrelation (Edson 2011)"
   Cuv = np.cov(signal.detrend(uvw_mc[:,0]),signal.detrend(uvwplat[:,0]))
   muuv = Cuv[0,1]/Cuv[1,1]
   uvw_mc[:,0] = uvw_mc[:,0] - muuv*signal.detrend(uvwplat[:,0])

   Cua = np.cov(signal.detrend(uvw_mc[:,0]),signal.detrend(acc[:,0]))
   muua = Cua[0,1]/Cua[1,1]
   uvw_mc[:,0] = uvw_mc[:,0] - muua*signal.detrend(acc[:,0])
   
   Cva = np.cov(signal.detrend(uvw_mc[:,1]),signal.detrend(acc[:,1]))
   muva = Cva[0,1]/Cva[1,1]
   uvw_mc[:,1] = uvw_mc[:,1] - muva*signal.detrend(acc[:,1])

   Cvv = np.cov(signal.detrend(uvw_mc[:,1]),signal.detrend(uvwplat[:,1]))
   muvv = Cvv[0,1]/Cvv[1,1]
   uvw_mc[:,1] = uvw_mc[:,1] - muvv*signal.detrend(uvwplat[:,1])
 
   Cwa = np.cov(signal.detrend(uvw_mc[:,2]),signal.detrend(acc[:,2]))
   muwa = Cwa[0,1]/Cwa[1,1]
   uvw_mc[:,2] = uvw_mc[:,2] - muwa*signal.detrend(acc[:,2])
   
   Cwv = np.cov(signal.detrend(uvw_mc[:,2]),signal.detrend(uvwplat[:,2]))
   muwv = Cwv[0,1]/Cwv[1,1]
   uvw_mc[:,2] = uvw_mc[:,2] - muwv*signal.detrend(uvwplat[:,2])
   #### END OF PART 4 ####

	#now reverse the sign of u to bring into meteorological convention
	uvw_mc[:,0] = -uvw_mc[:,0]

   #### PART 5 ####
   # Calculate fluxes, true wind speed, and averaging period statistics
   print "Calculating fluxes, true wind speed and associated statistics"

   # Calculate wind statistics and true wind speed & direction.
   # RelWdirGeo is apparent wind direction in earth coordinates
   #  (0 or 360 = N).  But this is still relative wind direction because
   # of effects of ship's velocity on the wind vector.
   U_mc = copy(uvw_mc[:,0])
   V_mc = copy(uvw_mc[:,1])
   W_mc = copy(uvw_mc[:,2])

   U_motcorr_ts = pd.Series(U_mc, index = timestamps)
   V_motcorr_ts = pd.Series(V_mc, index = timestamps)
   W_motcorr_ts = pd.Series(W_mc, index = timestamps)
   T_sonic_ts = pd.Series(T_K, index = timestamps)
   Zplat_ts = pd.Series(Zplat, index = timestamps)
   Wplat_ts = pd.Series(Wplat, index = timestamps)
   Aplat_ts = pd.Series(Aplat, index = timestamps)
   
   CO2_picarro_ts = pd.Series(CO2_picarro, index = timestamps)
   CH4_picarro_ts = pd.Series(CH4_picarro, index = timestamps)
   switch_state_picarro_ts = pd.Series(switch_state_picarro, index=timestamps)
   CO2_licor_ts = pd.Series(CO2_licor, index = timestamps)
   H2O_licor_ts = pd.Series(H2O_licor, index = timestamps)
   P_licor_ts = pd.Series(P_licor, index = timestamps)   
   sig_licor_ts = pd.Series(sig_licor, index = timestamps)
   
   heading_ts = pd.Series(heading, index = timestamps)
   SOS_ts = pd.Series(SOS, index = timestamps)
   ship_lat_ts = pd.Series(ship_lat, index = timestamps)
   ship_lon_ts = pd.Series(ship_lon, index = timestamps)

   U_motcorr_timestamps = U_motcorr_ts.resample('20min').mean()
   U_motcorr_timestamps = U_motcorr_timestamps.axes
   timestamp_limits = U_motcorr_timestamps[0]
   
   # Open file for writing output data
   if args.ofile:
      filename = args.ofile
      print "Output filename: ", filename
   else:
      filename = 'benchmark_20min.txt'
      
   if os.path.exists(filename):
      append_write = 'a' # append if already exists
      FILE = open(filename, append_write)
   else:
      append_write = 'w' # make a new file if not
      FILE = open(filename, append_write)
      # write a header
      FILE.write("StartDate StartTime EndDate EndTime Latitude Longitude ShipSpeed(knots) ShipBearing(deg) Utrue(m/s) RelWind(deg) Ustar(m/s) MomentumFlux(N/m2) SensibleHeatFlux(W/m2) LatentHeatFlux(W/m2) LongwaveFlux(W/m2) H2OConcMeanL(ppm) H2OConcStdL(ppm) SignalLmean(%) SignalLstd(%) CO2FluxP(mmol/m2/d) CO2ConcP(ppm) CO2ConcStdP(ppm) CO2FluxL(mmol/m2/d) CO2ConcL(ppm) CO2ConcStdL(ppm) CH4FluxP(umol/m2/d) CH4ConcP(ppm) CH4ConcStdP(ppm) AirTemp(degC) RelHum(%) TIR(W/m2) AtmosphericPressure(mb) SST(degC) Sal(PSU) RelWdirMean(deg) RelWdirStddev(deg) TrueWspdVar(m2/s2) Wvar(m2/s2) tilt(deg) UplatVar(m2/s2) VplatVar(m2/s2) WplatVar(m2/s2) HeadingRange(deg) HeadingSdev(deg) PicarroLag(s) PicarroResponse(s)\n") 

   for j in range(len(timestamp_limits) - 1):
      U_mc_sml = U_motcorr_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      V_mc_sml = V_motcorr_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      W_mc_sml = W_motcorr_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      T_sonic_sml = T_sonic_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      Zplat_sml = Zplat_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      Wplat_sml = Wplat_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      Aplat_sml = Aplat_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      
      CO2_picarro_sml = CO2_picarro_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      CH4_picarro_sml = CH4_picarro_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      switch_state_picarro_sml = switch_state_picarro_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      CO2_licor_sml = CO2_licor_ts[timestamp_limits[j]:timestamp_limits[j+1]] 
      H2O_licor_sml = H2O_licor_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      sig_licor_sml = sig_licor_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      P_licor_sml = P_licor_ts[timestamp_limits[j]:timestamp_limits[j+1]] 
      number_of_data_points = len(CO2_picarro_sml) # this is to check that enough data points in time frame
      
      if number_of_data_points < 1000:
         break
         
      heading_sml = heading_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      SOS_sml = SOS_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      Heading_Range = (max(heading_sml) - min(heading_sml))*rad2deg
      Heading_Sdev = np.std(heading_sml)*rad2deg
      
      lat_sml = ship_lat_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      lon_sml = ship_lon_ts[timestamp_limits[j]:timestamp_limits[j+1]]
      
      if args.underway:
         tair_sml = tair_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         rhum_sml = rhum_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         TIR_sml = TIR_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         atm_pr_sml = atm_pr_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         sst_sml = sst_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         sal_sml = sal_ts[timestamp_limits[j]:timestamp_limits[j+1]]
         tair = np.mean(tair_sml)
         rhum = np.mean(rhum_sml)
         TIR = np.mean(TIR_sml)
         atm_pr = np.mean(atm_pr_sml)
         sst = np.mean(sst_sml)
         sal = np.mean(sal_sml)

      UVW_mc_sml = np.zeros((len(U_mc_sml),3))*0.
      for k in range(len(U_mc_sml)):
         UVW_mc_sml[k,0] = U_mc_sml[k]
         UVW_mc_sml[k,1] = V_mc_sml[k]
         UVW_mc_sml[k,2] = W_mc_sml[k]
         
      # Correct winds for streamline distortion by ship superstructure.
      # Streamline corrected uvw (NOT corrected for ship speed) will be used
      # for flux calculations.
      # Function rotcorr.m rotates wind for zero hourly mean W and V.
      # Azimuth and Tilt are resultant rotation angles in degrees, and
      #    Azimuth should be equal to RelWdir after the rotation.
      # ** NOTE:  corrected uvw is no longer in ship coordinates. **
      # Vmean and Wmean should be zero after the rotation.
      uvw_str, tilt, azm, Umean, Vmean, Wmean = dbl_rot(UVW_mc_sml)

      #In vectors: Relative Wind = True Wind + ship velocity over water.
      #    Convert to rectangular coordinates, then 
      #    RelWind(X2,Y2) - SOS(X1,Y1) = TrueWind(X3,Y3)
      #    and then convert back to polar coordinates.
      # Utrue in this calculation is wind speed relative to the ocean
      #    surface. (i.e. accounts for ocean current)
      x1 = np.zeros(len(heading_sml))*0.
      x2 = np.zeros(len(heading_sml))*0.
      y1 = np.zeros(len(heading_sml))*0.
      y2 = np.zeros(len(heading_sml))*0.
      x3 = np.zeros(len(heading_sml))*0.
      y3 = np.zeros(len(heading_sml))*0.
      TrueWdir = np.zeros(len(heading_sml))*0.
      Utrue = np.zeros(len(heading_sml))*0.

      RelWdirGeo_sml = np.zeros(len(heading_sml))*0.
      RelWdir_sml = np.zeros(len(heading_sml))*0.
   
      for k in range(len(heading_sml)):
         #RelWdirGeo_sml[k] = rad2deg*np.arctan2(V_mc_sml[k],-1.*U_mc_sml[k])
         RelWdirGeo_sml[k] = rad2deg*np.arctan2(V_mc_sml[k],U_mc_sml[k])
         # output of atan2 in radians or degrees (zz = find(RelWdirGeo > 360))
         if (RelWdirGeo_sml[k] > 360.):
            RelWdirGeo_sml[k] = RelWdirGeo_sml[k] - 360.
         RelWdir_sml[k] = RelWdirGeo_sml[k] - (rad2deg*heading_sml[k]-360.)
         if (RelWdir_sml[k] > 180.):
            RelWdir_sml[k] = RelWdir_sml[k] - 360. # +/- 180 deg range 
      
         x1[k], y1[k] = pol2cart(rad2deg*heading_sml[k],SOS_sml[k])
         # uvw_str is streamline corrected uvw
         x2[k], y2[k] = pol2cart(RelWdirGeo_sml[k],uvw_str[k,0])
         x3[k] = x2[k] - x1[k]
         y3[k] = y2[k] - y1[k]
         TrueWdir[k], Utrue[k] = cart2pol(x3[k],y3[k]) # cart2pol returns direction in degrees
         if (TrueWdir[k] < 0): # check for negative angles
            TrueWdir[k] = TrueWdir[k] + 360.

      # Calculate statistics 
      RelWdirMean = np.mean(RelWdir_sml)
      RelWdirStddev = np.std(RelWdir_sml)
      TrueWdirMean = np.mean(TrueWdir)
      TrueWdirVar = np.var(TrueWdir)
      TrueWspdMean = np.mean(Utrue)
      TrueWspdVar = np.var(Utrue)
      ShipVelMean = np.mean(SOS_sml)
      ShipHeadingMean = np.mean(heading_sml)
      ShipLatMean = np.mean(lat_sml)
      ShipLonMean = np.mean(lon_sml)

      #### PART 6 ####
      # Use above information to work out EC fluxes 
      #   i. Momentum Flux
      #  ii. Sensible Heat Flux (SHF)
      # iii. Latent Heat Flux (LHF)
      #  iv. CO2 - Picarro
      #   v. CH4 - Picarro
      #  vi. CO2 - Licor
      # Calculate averaging period Ustar = sqrt(abs(mean(u'w')))
      Uprime = copy(uvw_str[:,0]) - Umean
      Wprime = copy(uvw_str[:,2]) - Wmean
      #Wprime_alt = alternative_Wprime(copy(uvw_str[:,2]))
      Ustar_cov = math.sqrt(math.fabs(np.mean(Uprime*Wprime)))
      # i. Momentum flux  [data source: sonic wind]
      rho_air = RHO_AIR
      momentum_flux = -1.0*rho_air*np.mean(Uprime*Wprime)
      
      # ii. Sensible Heat Flux [data source: sonic wind and sonic temperature]
      # Calculate the sensible heat flux using the (virtual) temperature from the sonic
      T_sonic_mean = np.mean(T_sonic_sml)
      Tprime = T_sonic_sml - T_sonic_mean
      SHF =  np.mean(Tprime*Wprime)*Cp_air*rho_air 

      # iii. Latent Heat Flux [data source: sonic wind and Licor H2O]
      H2O_licor_density = H2O_licor_sml # in mmol/m3

      if not(args.no_licor):
         dummy_licor_valve_switch_state = np.zeros(len(H2O_licor_sml))*0.
         q_licor_mean, q_licor_std, qF_licor = platmot_corrected_gas_flux(H2O_licor_density, dummy_licor_valve_switch_state, Zplat_sml, Wplat_sml, Aplat_sml, Wprime)
         # Units check for W/m2:
         # qF_licor: mmol/m2/s
         # Lvap:  J/g
         # M_w: g/mol      
         LHF = qF_licor*Lvap*M_w/mmol
      else:
         LHF = np.nan

      #  iv. CO2 Picarro
      #   v. CH4 Picarro
      #  vi. CO2 Licor
      # Units conversion
      # CO2 data from the Picarro and LICOR
      # CH4 data from Picarro
      #  convert from ppm to mmole/m3
      CO2_picarro_sml = mmol*CO2_picarro_sml*PRT/ppm
      CH4_picarro_sml = mmol*CH4_picarro_sml*PRT/ppm
         
      picarro_valve_switch_state = switch_state_picarro_sml
      CO2_picarro_mean, CO2_picarro_std, CO2F_picarro = platmot_corrected_gas_flux(CO2_picarro_sml, picarro_valve_switch_state, Zplat_sml, Wplat_sml, Aplat_sml, Wprime)
      CH4_picarro_mean, CH4_picarro_std, CH4F_picarro = platmot_corrected_gas_flux(CH4_picarro_sml, picarro_valve_switch_state, Zplat_sml, Wplat_sml, Aplat_sml, Wprime)

      #mmole/m2/d      
      CO2F_picarro = seconds_per_day*CO2F_picarro
      #umole/m2/d
      CH4F_picarro = umol*seconds_per_day*CH4F_picarro/mmol
      # back to ppm
      CO2_picarro_mean = CO2_picarro_mean*ppm/mmol/PRT
      CO2_picarro_std = CO2_picarro_std*ppm/mmol/PRT
      CH4_picarro_mean = CH4_picarro_mean*ppm/mmol/PRT
      CH4_picarro_std = CH4_picarro_std*ppm/mmol/PRT

      # Calculate the CO2 Flux using the data from the Licor
      if not(args.no_licor):
         sig_licor_mean = np.mean(sig_licor_sml)
         sig_licor_std = np.std(sig_licor_sml)

         # Calculate the mixing ratio of H2O in ppm based on a standard atmosphere
         H2O_licor_ppm = ppm*H2O_licor_density/mmol/PRT
         H2O_licor_mean = np.mean(H2O_licor_ppm)
         H2O_licor_std = np.std(H2O_licor_ppm)
         
         # calculate the flux of CO2 using the CO2 density measurements
         # convert to mmol/m3
         # CO2_licor_sml: mmol/m2/s
         dummy_licor_valve_switch_state = np.zeros(len(CO2_licor_sml))*0.
         CO2_licor_mean, CO2_licor_std, CO2F_licor = platmot_corrected_gas_flux(CO2_licor_sml, dummy_licor_valve_switch_state, Zplat_sml, Wplat_sml, Aplat_sml, Wprime)
         CO2F_licor = seconds_per_day*CO2F_licor
         
         # mixing ratio of CO2 in ppm needs to include the water vapour component
         CO2_licor_mean = 1e+3*CO2_licor_mean/(PRT - np.mean(H2O_licor_density)/1000.)
         CO2_licor_std = 1e+3*CO2_licor_std/(PRT - np.mean(H2O_licor_density)/1000.)
         # Calculate the longwave radiation flux (for completeness)
         P_licor_mean = np.mean(P_licor_sml)
         Fc = 0.50  # cloudiness correction factor: 1.0 is clear
         LWF = longwave_flux(T_sonic_mean,sst,Fc,P_licor_mean,H2O_licor_mean)
      else:
         H2O_licor_mean = np.nan 
         H2O_licor_std = np.nan
         sig_licor_mean = np.nan
         sig_licor_std = np.nan
         CO2_licor_mean = np.nan
         CO2_licor_std = np.nan
         CO2F_licor = np.nan
         LWF = np.nan
         
      # Calculate wind variances of rotated, motion corrected winds.
      # Again, these winds include ship velocity and are not in ship coord.
      Uvar = np.var(copy(uvw_str[:,0]))
      Vvar = np.var(copy(uvw_str[:,1]))
      Wvar = np.var(copy(uvw_str[:,2]))
      
      # Output of key parameters
      if (number_of_data_points > 1000):
         #print timestamp_limits[j], timestamp_limits[j+1], ShipLatMean, ShipLonMean, ShipVelMean/knots2ms, rad2deg*ShipHeadingMean, TrueWspdMean, RelWdirMean, \
         #Ustar_cov, momentum_flux, SHF, LHF, CO2F_picarro, CO2F_licor, CH4F_picarro, tair, rhum, TIR, atm_pr, sst, sal
      
         FILE.write("%s %s %7.3f %7.3f %7.3f %8.3f %7.3f %8.3f %8.3f %8.3f %8.3f %8.3f %8.3f %11.2f %8.3f %8.3f %8.3f %12.4e %8.3f %8.3f %12.4e %8.3f %8.3f %12.4e %11.6f %11.6f %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n" % (
         timestamp_limits[j], timestamp_limits[j+1], ShipLatMean, ShipLonMean, ShipVelMean/knots2ms, rad2deg*ShipHeadingMean, 
         TrueWspdMean, RelWdirMean, Ustar_cov, momentum_flux, SHF, LHF, LWF, H2O_licor_mean, H2O_licor_std, sig_licor_mean, sig_licor_std, CO2F_picarro, CO2_picarro_mean, CO2_picarro_std, 
         CO2F_licor, CO2_licor_mean, CO2_licor_std, CH4F_picarro, CH4_picarro_mean, CH4_picarro_std, tair, rhum, TIR, atm_pr, sst, sal,
         math.fabs(RelWdirMean), RelWdirStddev, TrueWspdVar, Wvar, tilt, uplat_Var, vplat_Var, wplat_Var,
         Heading_Range, Heading_Sdev, Picarro_lag_s, Picarro_response_s))
      
   FILE.close()

   sys.exit()
      
if __name__=='__main__':
   main()
   
