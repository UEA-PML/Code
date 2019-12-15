import sys
import datetime as dt
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import dates
from matplotlib import rc
from datetime import datetime
from datetime import date, timedelta
from numpy import copy, nanmean
import pylab 
import math
import argparse
import pandas as pd
import time

# Function: convert_as_needed()
def convert_as_needed(ts):
    dat_time = ''
    sec_frac = ts[0]
    try:
        # parse strings
        dat_time = datetime.strptime(' '.join(ts), '%Y-%m-%d %H:%M:%S.%f')
    except:
        # assume doesn't have decimal fraction of a second
        sec_frac = ts[0]
        sec_frac = sec_frac + ".000001"
        dat_time = datetime.strptime(sec_frac, '%Y-%m-%d %H:%M:%S.%f')

    return dat_time

# Function: read_underway(aco)
# Author: Tim Smyth
# Description:
#  Function to read in JCR underway data
# Inputs:
#  ASCII File - formatted and cleaned using pre-processing shell script
#             - /tmp/process_ecco2flux_oceanlogger.txt
# Outputs:
#  1) Floating point matrices: 
#  tair - air temperature (degC)
#  rhum - relative humidity (%)
#  TIR - Total Irradiance (W/m2)
#  atm_pr - atmospheric pressure (mb)
#  sst - sea surface temperature (degC)
#  sal - salinity (PSU)
# 
#  2) List of timestamps: datetimeList
# 
# Call function as:
#  tair, rhum, TIR, atm_pr, sst, sal, datetimeList = read_underway(aco)
#
def read_underway(aco):
   with open("/tmp/process_ecco2flux_oceanlogger.txt") as f:
      underway_data = f.readlines()
   
   tt = [row.split(' ')[0:2] for row in underway_data]

   if aco:   
      yyyy = [row.split(' ')[0] for row in underway_data]
      day_frac = [row.split(' ')[1] for row in underway_data]
   
      fday_frac = np.zeros(len(day_frac), dtype=float)

      for dd in range(len(day_frac)):
         fday_frac[dd] = float(day_frac[dd])
         frac = math.modf(fday_frac[dd])
         hour = math.modf(frac[0]*24)
         minute = math.modf(hour[0]*60)
         second = math.modf(minute[0]*60)
         tt[dd] = [yyyy[dd], str(int(frac[1])), str(int(hour[1])), str(int(minute[1])), str(int(second[1]))]
     
      datetimeList = [datetime.strptime(' '.join(row),'%Y %j %H %M %S') for row in tt]

   else: 
      datetimeList = [datetime.strptime(' '.join(row),'%m/%d/%Y %H:%M:%S.%f') for row in tt]


   tair1 = [row.split(' ')[2] for row in underway_data]
   rhum1 = [row.split(' ')[3] for row in underway_data]
   par1 = [row.split(' ')[4] for row in underway_data]
   TIR1 = [row.split(' ')[5] for row in underway_data]
   tair2 = [row.split(' ')[6] for row in underway_data]
   rhum2 = [row.split(' ')[7] for row in underway_data]
   par2 = [row.split(' ')[8] for row in underway_data]
   TIR2 = [row.split(' ')[9] for row in underway_data]
   atm_pr1 = [row.split(' ')[10] for row in underway_data]
   atm_pr2 = [row.split(' ')[11] for row in underway_data]
   seatemp = [row.split(' ')[12] for row in underway_data]
   salinity = [row.split(' ')[13] for row in underway_data]
  
   tair = np.zeros(len(tair1))*0.
   rhum = np.zeros(len(rhum1))*0.
   TIR = np.zeros(len(par1))*0.
   atm_pr = np.zeros(len(atm_pr1))*0.
   sst = np.zeros(len(seatemp))*0.
   sal = np.zeros(len(salinity))*0.

   for j in range(len(tair1)):
      # Ensure temperature and humidity as shaded as possible
      tair[j] = np.min([float(tair1[j]),float(tair2[j])])
      rhum[j] = np.min([float(rhum1[j]),float(rhum2[j])])
      # Ensure non-shading of TIR
      TIR[j] = np.max([float(TIR1[j]),float(TIR2[j])])
      # Set very low values of TIR to zero
      if (TIR[j] < 1.):
         TIR[j] = 0.
      atm_pr[j] = np.mean([float(atm_pr1[j]),float(atm_pr2[j])])
      sst[j] = float(seatemp[j])
      sal[j] = float(salinity[j])

   return tair, rhum, TIR, atm_pr, sst, sal, datetimeList      


# Function: read_licor()
# Author: Tim Smyth
# Description:
#  Function to read in Licor data
# Inputs:
#  ASCII File - formatted and cleaned using pre-processing shell script
#             - /tmp/process_ecco2flux_licor.txt
# Outputs:
#  1) Floating point matrices: fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag
#   fco2_density: N element array of CO2 density
#   fh2o_density: N element array of H2O density
#   fco2_ppm: N element array of CO2 concentration in ppm
#   fh2o_ppm: N element array of H2O concentration in ppm
#   ftemp: N element array of Licor temperature in degC
#   fpressure: N element array of atmospheric pressure in hPa
#   fdiag: N element array of diagnosic variable
#  2) List of timestamps: datetimeList
# 
# Call function as:
#  fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList_licor = read_licor()
#
def read_licor():
   with open("/tmp/process_ecco2flux_licor.txt") as f:
      licor_data = f.readlines()
   
   tt = [row.split(' ')[0:5] for row in licor_data]
   diag = [row.split(' ')[9] for row in licor_data]
   co2_raw = [row.split(' ')[11] for row in licor_data]
   co2_density = [row.split(' ')[13] for row in licor_data]
   h2o_raw = [row.split(' ')[15] for row in licor_data]
   h2o_density = [row.split(' ')[17] for row in licor_data]
   temp = [row.split(' ')[19] for row in licor_data]
   pressure = [row.split(' ')[21] for row in licor_data]
   
   datetimeList = [datetime.strptime(' '.join(row),'%a %b %d %H:%M:%S.%f %Y') for row in tt]
   
   idiag = np.zeros(len(diag), dtype=int)
   fdiag = np.zeros(len(diag))*0.
   fco2_raw = np.zeros(len(co2_raw))*0.
   fh2o_raw = np.zeros(len(h2o_raw))*0.
   fco2_density = np.zeros(len(co2_density))*0.
   fh2o_density = np.zeros(len(h2o_density))*0.
   ftemp = np.zeros(len(temp))*0.
   fpressure = np.zeros(len(pressure))*0.
   fPRT = np.zeros(len(pressure))*0.
   fco2_ppm = np.zeros(len(co2_density))*0.
   fh2o_ppm = np.zeros(len(h2o_density))*0.
   
   for j in range(len(diag)):
      idiag[j] = int(diag[j])
      # Convert decimal (integer) into binary
      bdiag = "{0:07b}".format(idiag[j])
      fdiag[j] = int(bdiag[4:8],2)*6.67
      fco2_raw[j] = float(co2_raw[j])
      fh2o_raw[j] = float(h2o_raw[j])
      fco2_density[j] = float(co2_density[j])
      fh2o_density[j] = float(h2o_density[j])
      ftemp[j] = float(temp[j])
      fPRT[j] = (float(pressure[j])*1000.)/(8.31451*(ftemp[j]+273.15))
      fpressure[j] = float(pressure[j])
      fco2_ppm[j] = 1e+3*fco2_density[j]/(fPRT[j] - fh2o_density[j]/1000.)
      fh2o_ppm[j] = 1e+3*fh2o_density[j]/fPRT[j]
   
   return fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList

# Function: read_seatex_gga()
# Author: Tim Smyth
# Description: 
#  Read in a seatex gga data file which contains timestamped latitude and longitude 
#
# Inputs:
#  ASCII file - formatted and cleaned by pre-processing shell script
#             - /tmp/process_ecco2flux_seatex-gga.txt
# Outputs:
#  latitude_dec - decimal latitude, (-ve S, +ve N)
#  longitude_dec - decimal longitude (-ve W, +ve E) 
#  datetimeList - timestamp
#
# Call function as:
#  lat, lon, datetimeList_seatex_gga = read_seatex_gga()
def read_seatex_gga():
   with open("/tmp/process_ecco2flux_seatex-gga.txt") as f:
      location_data = f.readlines()
   
   tt = [row.split(' ')[0:2] for row in location_data]
   yyyy = [row.split(' ')[0] for row in location_data]
   day_frac = [row.split(' ')[1] for row in location_data]
   latitude = [row.split(' ')[2] for row in location_data]
   longitude = [row.split(' ')[3] for row in location_data]
   
   fday_frac = np.zeros(len(day_frac), dtype=float)

   for dd in range(len(day_frac)):
     fday_frac[dd] = float(day_frac[dd])
     frac = math.modf(fday_frac[dd])
     hour = math.modf(frac[0]*24)
     minute = math.modf(hour[0]*60)
     second = math.modf(minute[0]*60)
     tt[dd] = [yyyy[dd], str(int(frac[1])), str(int(hour[1])), str(int(minute[1])), str(int(second[1]))]
     
   datetimeList = [datetime.strptime(' '.join(row),'%Y %j %H %M %S') for row in tt]
   
   latitude_dec = np.zeros(len(latitude))*0.
   longitude_dec = np.zeros(len(longitude))*0.
   
   for j in range(len(latitude)):
      latitude_dec[j] = float(latitude[j]) 
      longitude_dec[j] = float(longitude[j])

   return latitude_dec, longitude_dec, datetimeList

# Function: read_seatex_gll()
# Author: Tim Smyth
# Description: 
#  Read in a seatex gll data file which contains timestamped latitude and longitude in NMEA format
#
# Inputs:
#  ASCII file - formatted and cleaned by pre-processing shell script
#             - /tmp/process_ecco2flux_seatex-gll.txt
# Outputs:
#  latitude_dec - decimal latitude, (-ve S, +ve N)
#  longitude_dec - decimal longitude (-ve W, +ve E) 
#  datetimeList - timestamp
#
# Call function as:
#  lat, lon, datetimeList_seatex_gll = read_seatex_gll()
def read_seatex_gll():
   with open("/tmp/process_ecco2flux_seatex-gll.txt") as f:
      location_data = f.readlines()
   
   tt = [row.split(' ')[0:2] for row in location_data]
   latitude = [row.split(' ')[2] for row in location_data]
   latitude_NS = [row.split(' ')[3] for row in location_data]
   longitude = [row.split(' ')[4] for row in location_data]
   longitude_EW = [row.split(' ')[5] for row in location_data]

   datetimeList = [datetime.strptime(' '.join(row),'%m/%d/%Y %H:%M:%S.%f') for row in tt]
   
   latitude_dec = np.zeros(len(latitude))*0.
   longitude_dec = np.zeros(len(longitude))*0.
   
   for j in range(len(latitude)):
      # convert from NMEA format DDMM.MMMMM to DD.ddddd
      latitude_deg = int(float(latitude[j])/100.)
      latitude_minutes = float(latitude[j]) - latitude_deg*100.
      latitude_dec[j] = latitude_deg + latitude_minutes/60.
      if (latitude_NS[j] == 'S'):
         latitude_dec[j] = -1.*latitude_dec[j]
      
      longitude_deg = int(float(longitude[j])/100.)
      longitude_minutes = float(longitude[j]) - longitude_deg*100.
      longitude_dec[j] = longitude_deg + longitude_minutes/60.
      #print longitude_EW[j]
      if (longitude_EW[j] == 'W\n'): # Take care here as this includes a carriage return which has not been eradicated
         longitude_dec[j] = -1.*longitude_dec[j]
      
   return latitude_dec, longitude_dec, datetimeList

# Function: read_seatex_vtg(aco)
# Author: Tim Smyth
# Description: 
#  Read in a seatex vtg data file which contains timestamped ships heading and speed
#
# Inputs:
#  ASCII file - formatted and cleaned by pre-processing shell script
#             - /tmp/process_ecco2flux_seatex-vtg.txt
# Outputs:
#  fship_heading - N element array of floating point ship heading (degrees)
#  fship_speed - N element array of floating point ship speed (knots)
#  datetimeList - timestamp
#
# Call function as:
#  heading, SOG, datetimeList_seatex_vtg = read_seatex_vtg(aco)
def read_seatex_vtg(aco):
   with open("/tmp/process_ecco2flux_seatex-vtg.txt") as f:
      ship_nav_data = f.readlines()
     
   if aco:   
      tt = [row.split(' ')[0:2] for row in ship_nav_data]
      yyyy = [row.split(' ')[0] for row in ship_nav_data]
      day_frac = [row.split(' ')[1] for row in ship_nav_data]
      ship_heading = [row.split(' ')[2] for row in ship_nav_data]
      ship_speed = [row.split(' ')[3] for row in ship_nav_data]
   
      fday_frac = np.zeros(len(day_frac), dtype=float)

      for dd in range(len(day_frac)):
         fday_frac[dd] = float(day_frac[dd])
         frac = math.modf(fday_frac[dd])
         hour = math.modf(frac[0]*24)
         minute = math.modf(hour[0]*60)
         second = math.modf(minute[0]*60)
         tt[dd] = [yyyy[dd], str(int(frac[1])), str(int(hour[1])), str(int(minute[1])), str(int(second[1]))]
     
      datetimeList = [datetime.strptime(' '.join(row),'%Y %j %H %M %S') for row in tt]
   else: 
      tt = [row.split(' ')[0:2] for row in ship_nav_data]
      ship_heading = [row.split(' ')[2] for row in ship_nav_data]
      ship_speed = [row.split(' ')[3] for row in ship_nav_data]
      datetimeList = [datetime.strptime(' '.join(row),'%m/%d/%Y %H:%M:%S.%f') for row in tt]

   fship_heading = np.zeros(len(ship_heading))*0.
   fship_speed = np.zeros(len(ship_speed))*0.
   
   for j in range(len(ship_heading)):
      fship_heading[j] = float(ship_heading[j])
      fship_speed[j] = float(ship_speed[j])
   
   return fship_heading, fship_speed, datetimeList

# Function: read_gyro(aco)
# Author: Tim Smyth
# Description: 
#  Read in a gyro data file which contains timestamped ships heading 
#
# Inputs:
#  ASCII file - formatted and cleaned by pre-processing shell script
#             - /tmp/process_ecco2flux_gyro.txt
# Outputs:
#  fship_heading - N element array of floating point ship heading (degrees)
#  datetimeList - timestamp
#
# Call function as:
#  heading, datetimeList_gyro = read_gyro()
def read_gyro(aco):
   with open("/tmp/process_ecco2flux_gyro.txt") as f:
      ship_gyro_data = f.readlines()
      
   tt = [row.split(' ')[0:2] for row in ship_gyro_data]

   if aco:
      yyyy = [row.split(' ')[0] for row in ship_gyro_data]
      day_frac = [row.split(' ')[1] for row in ship_gyro_data]
      ship_heading = [row.split(' ')[2] for row in ship_gyro_data]
   
      fday_frac = np.zeros(len(day_frac), dtype=float)

      for dd in range(len(day_frac)):
         fday_frac[dd] = float(day_frac[dd])
         frac = math.modf(fday_frac[dd])
         hour = math.modf(frac[0]*24)
         minute = math.modf(hour[0]*60)
         second = math.modf(minute[0]*60)
         tt[dd] = [yyyy[dd], str(int(frac[1])), str(int(hour[1])), str(int(minute[1])), str(int(second[1]))]
     
      datetimeList = [datetime.strptime(' '.join(row),'%Y %j %H %M %S') for row in tt]

   else: 
      ship_heading = [row.split(' ')[2] for row in ship_gyro_data]
      datetimeList = [datetime.strptime(' '.join(row),'%m/%d/%Y %H:%M:%S.%f') for row in tt]

   fship_heading = np.zeros(len(ship_heading))*0.
   
   for j in range(len(ship_heading)):
      fship_heading[j] = float(ship_heading[j])
   
   return fship_heading, datetimeList
   
# Function: read_LPMS()
# Author: Tim Smyth
# Description: 
#  Function to read in LPMS data
#
# Inputs:
#   ASCII file - cleaned and reformatted using pre-processing shell script
#              - /tmp/process_ecco2flux_LPMS.txt
# Outputs:
#  frotx_degs - N element array of x-axis rotation in degrees 
#  froty_degs - N element array of y-axis rotation in degrees 
#  frotz_degs - N element array of z-axis rotation in degrees 
#  faccelx_g  - N element array of x-axis acceleration in 'g' (i.e. 1g = 9.8 m/s2) 
#  faccely_g  - N element array of y-axis acceleration in 'g' 
#  faccelz_g  - N element array of z-axis acceleration in 'g'
#  DatetimeList - list of timestamps
#
# Call function as:
#  lrotx_degs, lroty_degs, lrotz_degs, laccelx_g, laccely_g, laccelz_g, datetimeList_LPMS = read_LPMS()
def read_LPMS():
   with open("/tmp/process_ecco2flux_LPMS.txt") as f:
      motion_data = f.readlines()

   deg2rad = math.pi/180.
   rad2deg = 180./math.pi
   
   tt = [row.split(',')[0:6] for row in motion_data]
   seconds = [row.split(',')[6] for row in motion_data]
   rotx = [row.split(',')[7] for row in motion_data]
   roty = [row.split(',')[8] for row in motion_data]
   rotz = [row.split(',')[9] for row in motion_data]
   accelx = [row.split(',')[10] for row in motion_data] 
   accely = [row.split(',')[11] for row in motion_data] 
   accelz = [row.split(',')[12] for row in motion_data]
   
   datetimeList = [datetime.strptime(' '.join(row),'%Y %m %d %H %M %S.%f') for row in tt]

   frotx_degs = np.zeros(len(rotx))*0.
   froty_degs = np.zeros(len(roty))*0.
   frotz_degs = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_degs[j] = rad2deg*float(rotx[j])/1000.
      froty_degs[j] = rad2deg*float(roty[j])/1000.
      frotz_degs[j] = rad2deg*float(rotz[j])/1000.
      faccelx_g[j] = float(accelx[j])/-1000.
      faccely_g[j] = float(accely[j])/-1000.
      faccelz_g[j] = float(accelz[j])/-1000.      

   return frotx_degs, froty_degs, frotz_degs, faccelx_g, faccely_g, faccelz_g, datetimeList

# Function: read_CR6()
# Author: Tim Smyth
# Description: 
#  Function to read in data from CR6 logger
#
# Inputs:
#   ASCII file - cleaned and reformatted using pre-processing shell script
#              - /tmp/process_ecco2flux_CR6.txt
# Outputs:
#  frotx_degs - N element array of x-axis rotation in degrees 
#  froty_degs - N element array of y-axis rotation in degrees 
#  frotz_degs - N element array of z-axis rotation in degrees 
#  faccelx_g  - N element array of x-axis acceleration in 'g' (i.e. 1g = 9.8 m/s2) 
#  faccely_g  - N element array of y-axis acceleration in 'g' 
#  faccelz_g  - N element array of z-axis acceleration in 'g'
#  u_ms - N element array of floating point horizontal vector wind (u) in m/s
#  v_ms - N element array of floating point horizontal vector wind (v) in m/s
#  w_ms - N element array of floating point vertical vector wind (w) in m/s
#  fco2_density - N element array of CO2 density
#  fh2o_density - N element array of H2O density
#  fco2_ppm - N element array of CO2 concentration in ppm
#  fh2o_ppm - N element array of H2O concentration in ppm
#  ftemp - N element array of Licor temperature in degC
#  fpressure - N element array of atmospheric pressure in hPa
#  fdiag - N element array of diagnosic variable
#  datetimeList - list of timestamps
#
# Call function as:
#   u_ms, v_ms, w_ms, t_degC, lrotx_degs, lroty_degs, lrotz_degs, laccelx_g, laccely_g, laccelz_g, fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList_CR6 = read_CR6()
def read_CR6():
   with open("/tmp/process_ecco2flux_CR6.txt") as f:
      CR6_data = f.readlines()

   deg2rad = math.pi/180.
   rad2deg = 180./math.pi

   # Instrument order is: sonic, LPMS, Licor
   # Specifics are: "TIMESTAMP","SonicX","SonicY","SonicZ","SonicT",\
   #                "RotX","RotY","RotZ","AccX","AccY","AccZ",\
   #                "LicorCO2D","LicorH2OD","LicorTemp","LicorPress","LicorDiag"

   tt = [row.split(',')[0:1] for row in CR6_data]
   
   # wind data
   x = [row.split(',')[1] for row in CR6_data]
   y = [row.split(',')[2] for row in CR6_data]
   z = [row.split(',')[3] for row in CR6_data]
   t = [row.split(',')[4] for row in CR6_data]

   u_ms = np.zeros(len(y))*0.
   v_ms = np.zeros(len(x))*0.
   w_ms = np.zeros(len(z))*0.
   t_degC = np.zeros(len(t))*0.
   
   for j in range(len(x)):
      if x[j] != '':
         # Change of order x -> -v; y-> u
         #u_ms[j] = float(x[j])/100.
         #v_ms[j] = float(y[j])/100.
         v_ms[j] = float(x[j])/-100.
         u_ms[j] = float(y[j])/100.
         w_ms[j] = float(z[j])/100.
         t_degC[j] = float(t[j])/100.
      else:
         v_ms[j] = float('NaN')
         u_ms[j] = float('NaN')
         w_ms[j] = float('NaN')
         t_degC[j] = float('NaN')

   # motion data
   rotx = [row.split(',')[5] for row in CR6_data]
   roty = [row.split(',')[6] for row in CR6_data]
   rotz = [row.split(',')[7] for row in CR6_data]
   accelx = [row.split(',')[8] for row in CR6_data] 
   accely = [row.split(',')[9] for row in CR6_data] 
   accelz = [row.split(',')[10] for row in CR6_data]

   frotx_degs = np.zeros(len(rotx))*0.
   froty_degs = np.zeros(len(roty))*0.
   frotz_degs = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_degs[j] = rad2deg*float(rotx[j])/1000.
      froty_degs[j] = rad2deg*float(roty[j])/1000.
      frotz_degs[j] = rad2deg*float(rotz[j])/1000.
      faccelx_g[j] = float(accelx[j])/-1000.
      faccely_g[j] = float(accely[j])/-1000.
      faccelz_g[j] = float(accelz[j])/-1000.      

   # Licor data
   co2_density = [row.split(',')[11] for row in CR6_data]
   h2o_density = [row.split(',')[12] for row in CR6_data]
   temp = [row.split(',')[13] for row in CR6_data]
   pressure = [row.split(',')[14] for row in CR6_data]
   diag = [row.split(',')[15] for row in CR6_data]
   
   #datetimeList = [datetime.strptime(' '.join(row),'%Y-%m-%d %H:%M:%S.%f') for row in tt]
   datetimeList = [convert_as_needed(row) for row in tt]

   idiag = np.zeros(len(diag), dtype=int)
   fdiag = np.zeros(len(diag))*0.
   fco2_density = np.zeros(len(co2_density))*0.
   fh2o_density = np.zeros(len(h2o_density))*0.
   ftemp = np.zeros(len(temp))*0.
   fpressure = np.zeros(len(pressure))*0.
   fPRT = np.zeros(len(pressure))*0.
   fco2_ppm = np.zeros(len(co2_density))*0.
   fh2o_ppm = np.zeros(len(h2o_density))*0.

   for j in range(len(diag)):
      diag_no_cr = diag[j].rstrip()
      if (diag_no_cr == "NAN"):
         diag_no_cr = "0"
      idiag[j] = int(float(diag_no_cr))
      # Convert decimal (integer) into binary
      bdiag = "{0:07b}".format(idiag[j])
      fdiag[j] = int(bdiag[4:8],2)*6.67
      fco2_density[j] = float(co2_density[j])
      fh2o_density[j] = float(h2o_density[j])
      ftemp[j] = float(temp[j])
      fPRT[j] = (float(pressure[j])*1000.)/(8.31451*(ftemp[j]+273.15))
      fpressure[j] = float(pressure[j])*10.#/101.325
      fco2_ppm[j] = 1e+3*fco2_density[j]/(fPRT[j] - fh2o_density[j]/1000.)
      fh2o_ppm[j] = 1e+3*fh2o_density[j]/fPRT[j]
      
      if (fdiag[j] <= 40.0):
         fco2_ppm[j] = np.nan
         fh2o_ppm[j] = np.nan
         fdiag[j] = np.nan
         ftemp[j] = np.nan
         fpressure[j] = np.nan

   return u_ms, v_ms, w_ms, t_degC, frotx_degs, froty_degs, frotz_degs, faccelx_g, faccely_g, faccelz_g, fco2_density, fh2o_density, fco2_ppm, fh2o_ppm, ftemp, fpressure, fdiag, datetimeList

# Function: read_picarro()
# Author: Tim Smyth
# Description: 
#  Function to read in Picarro data
#
# Inputs:
#  ASCII file - cleaned and formatted using pre-processing shell script
#             - /tmp/process_ecco2flux_picarro.txt
# Outputs:
#  fch4 - N element array of floating point CH4 cocentration (ppm)
#  fco2_dry - N element array of floating point CO2 concentration (ppm)
#  fswitch_state - N element array of switch state
#  fcavity_pressure_torr - N element array of floating point cavity pressure (atmospheres)
#  datetimeList - List of timestamps
#
# Call function as:
#  fch4, fco2_dry, fswitch_state, fcavity_pressure_atm_co2, datetimeList = read_picarro()
def read_picarro():
   with open("/tmp/process_ecco2flux_picarro.txt") as f:
      picarro_data = f.readlines()
   
   tt = [row.split(' ')[0:5] for row in picarro_data]
   seconds = [row.split(' ')[5] for row in picarro_data]
   cavity_pressure_torr = [row.split(' ')[6] for row in picarro_data]
   co2_temperature = [row.split(' ')[7] for row in picarro_data]
   switch_state = [row.split(' ')[8] for row in picarro_data]
   
   ch4 = [row.split(' ')[9] for row in picarro_data]
   co2_air = [row.split(' ')[10] for row in picarro_data]
   co2_dry = [row.split(' ')[11] for row in picarro_data]
   h2o = [row.split(' ')[12] for row in picarro_data]

   # datetimeList based on the timestamp from the PC
   datetimeListPC = [datetime.strptime(' '.join(row),'%a %b %d %H:%M:%S.%f %Y') for row in tt]
   
   # Alternative datetimeList based on Picarro timestamp
   tts = []
   for j in range(len(seconds)):
       mlsec = repr(float(seconds[j])).split('.')[1][:3]
       tts.append(time.strftime('%a %b %d %H:%M:%S.{} %Y'.format(mlsec), time.localtime(float(seconds[j]))))    
   datetimeList = [datetime.strptime(row,'%a %b %d %H:%M:%S.%f %Y') for row in tts]
   
   fcavity_pressure_torr = np.zeros(len(cavity_pressure_torr))*0.
   fco2_temperature = np.zeros(len(co2_temperature))*0.
   fch4 = np.zeros(len(ch4))*0.
   fco2_air = np.zeros(len(co2_air))*0.
   fco2_dry = np.zeros(len(co2_dry))*0.
   fh2o = np.zeros(len(h2o))*0.
   fswitch_state = np.zeros(len(switch_state))*0.
   
   for j in range(len(ch4)):
      fch4[j] = float(ch4[j])
      fco2_air[j] = float(co2_air[j])
      fco2_dry[j] = float(co2_dry[j])
      fh2o[j] = float(h2o[j])
      fcavity_pressure_torr[j] = float(cavity_pressure_torr[j])
      fco2_temperature[j] = float(co2_temperature[j])
      fswitch_state[j] = float(switch_state[j]) - 32.0

   # Work out the decorrelation between the cavity pressure and the gas (CO2 and CH4)
   # rather than doing a simple mean, produce a linear fit
   # Create an array of same size as input_gas which increments at 0.1
   time_s_array = np.zeros(len(fco2_air))*0.
   for ss in range(len(fco2_air)):
      time_s_array[ss] = float(ss)*0.1

   m_press,press_c = np.polyfit(time_s_array,fcavity_pressure_torr,1)
   press_lin = press_c + m_press*time_s_array
   press_prime = fcavity_pressure_torr - press_lin

   m_co2,co2_c = np.polyfit(time_s_array,fco2_dry,1)
   co2_lin = co2_c + m_co2*time_s_array
   co2_prime = fco2_dry - co2_lin

   m_ch4,ch4_c = np.polyfit(time_s_array,fch4,1)
   ch4_lin = ch4_c + m_ch4*time_s_array
   ch4_prime = fch4 - ch4_lin
   
   press_cov = co2_prime*press_prime
   mu_CO2 = np.mean(press_cov)/np.var(press_prime)

   press_cov = ch4_prime*press_prime
   mu_ch4 = np.mean(press_cov)/np.var(press_prime)

   fco2_dry = fco2_dry - mu_CO2*press_prime
   fch4 = fch4 - mu_ch4*press_prime
   
   return fch4, fco2_dry, fswitch_state, fcavity_pressure_torr, datetimeList
   
# Function: read_sysdon()
# Author: Tim Smyth
# Description: 
#  Function to read in SysDon data
#
# Inputs:
#  ASCII file - cleaned and formatted using pre-processing shell script
#             - /tmp/process_ecco2flux_SysDon.txt
# Outputs:
#  frotx_degs   - N element array of x-axis rotation in degrees 
#  froty_degs   - N element array of y-axis rotation in degrees 
#  frotz_degs   - N element array of z-axis rotation in degrees 
#  faccelx_g    - N element array of x-axis acceleration in 'g' (i.e. 1g = 9.8 m/s2) 
#  faccely_g    - N element array of y-axis acceleration in 'g' 
#  faccelz_g    - N element array of z-axis acceleration in 'g'
#  datetimeList - list of timestamps
#
# Call function as:
#  frotx_degs, froty_degs, frotz_degs, faccelx_g, faccely_g, faccelz_g, datetimeList = read_sysdon()
def read_sysdon():
   with open("/tmp/process_ecco2flux_SysDon.txt") as f:
      motion_data = f.readlines()
   
   seconds = [row.split('\t')[0] for row in motion_data]
   accelz = [row.split('\t')[2] for row in motion_data] 
   accely = [row.split('\t')[3] for row in motion_data] 
   accelx = [row.split('\t')[4] for row in motion_data]
   rotz = [row.split('\t')[5] for row in motion_data]
   roty = [row.split('\t')[6] for row in motion_data]
   rotx = [row.split('\t')[7] for row in motion_data]

   tt = []
   
   date0 = datetime(1904,1,1)
   for j in range(len(seconds)):
      delta = timedelta(seconds=float(seconds[j])+0.000001)
      tt.append(str(date0 + delta))

   datetimeList = [datetime.strptime(''.join(row),'%Y-%m-%d %H:%M:%S.%f') for row in tt]   

   frotx_degs = np.zeros(len(rotx))*0.
   froty_degs = np.zeros(len(roty))*0.
   frotz_degs = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_degs[j] = -1.0*(((float(rotx[j])*5./4095.)-2.5)/0.027)
      froty_degs[j] = -1.0*(((float(roty[j])*5./4095.)-2.5)/0.027)
      frotz_degs[j] = -1.0*(((float(rotz[j])*5./4095.)-2.5)/0.027)
      faccelx_g[j] = 1.0*(((float(accelx[j])*5./4095.)-2.5)/0.75)
      faccely_g[j] = 1.0*(((float(accely[j])*5./4095.)-2.5)/0.75)
      faccelz_g[j] = 1.0*(((float(accelz[j])*5./4095.)-2.5)/0.75)      

   return frotx_degs, froty_degs, frotz_degs, faccelx_g, faccely_g, faccelz_g, datetimeList

# Function: read_wind()
# Author: Tim Smyth
# Description: 
#  Function to read in wind data from Sonic 3D-wind anemometer
#
# Inputs:
#  ASCII file - cleaned and formatted using pre-processing shell script
#             - /tmp/process_ecco2flux_wind.txt
# Outputs:
#  u_ms - N element array of floating point horizontal vector wind (u) in m/s
#  v_ms - N element array of floating point horizontal vector wind (v) in m/s
#  w_ms - N element array of floating point vertical vector wind (w) in m/s
#  t_degC - N element array of floating point temperature (degC)
#  datetimeList - list of timestamps
#
# Call function as:
# u_ms, v_ms, w_ms, t_degC, datetimeList = read_wind()    
def read_wind():
   with open("/tmp/process_ecco2flux_wind.txt") as f:
      wind_data = f.readlines()
   
   tt = [row.split(' ')[0:5] for row in wind_data]
   x  = [row.split(' ')[5] for row in wind_data]
   y  = [row.split(' ')[6] for row in wind_data]
   z  = [row.split(' ')[7] for row in wind_data]
   t  = [row.split(' ')[8] for row in wind_data]
   
   datetimeList = [datetime.strptime(' '.join(row),'%a %b %d %H:%M:%S.%f %Y') for row in tt]
   
   u_ms = np.zeros(len(y))*0.
   v_ms = np.zeros(len(x))*0.
   w_ms = np.zeros(len(z))*0.
   t_degC = np.zeros(len(t))*0.
   
   for j in range(len(x)):
      if x[j] != '':
         # Change of order x -> -v; y-> u
         #u_ms[j] = float(x[j])/100.
         #v_ms[j] = float(y[j])/100.
         v_ms[j] = float(x[j])/-100.
         u_ms[j] = float(y[j])/100.
         w_ms[j] = float(z[j])/100.
         t_degC[j] = float(t[j])/100.
      else:
         v_ms[j] = float('NaN')
         u_ms[j] = float('NaN')
         w_ms[j] = float('NaN')
         t_degC[j] = float('NaN')

   return u_ms, v_ms, w_ms, t_degC, datetimeList

# Function: read_gill()
# Author: Tim Smyth
# Description: 
#  Function to read in wind data from Gill 3D-wind anemometer
#
# Inputs:
#  ASCII file - cleaned and formatted using pre-processing shell script
#             - /tmp/process_ecco2flux_gill.txt
# Outputs:
#  u_ms - N element array of floating point horizontal vector wind (u) in m/s
#  v_ms - N element array of floating point horizontal vector wind (v) in m/s
#  w_ms - N element array of floating point vertical vector wind (w) in m/s
#  datetimeList - list of timestamps
#
# Call function as:
# u_ms, v_ms, w_ms, datetimeList = read_gill()    
def read_gill():
   with open("/tmp/process_ecco2flux_gill.txt") as f:
      wind_data = f.readlines()
      
   tt = [row.split(',')[0:1] for row in wind_data]
   u = [row.split(',')[2] for row in wind_data]
   v = [row.split(',')[3] for row in wind_data]
   w = [row.split(',')[4] for row in wind_data]

   datetimeList = [datetime.strptime(','.join(row),'%Y%m%dT%H%M%SZ') for row in tt]
   
   u_ms = np.zeros(len(u))*0.
   v_ms = np.zeros(len(v))*0.
   w_ms = np.zeros(len(w))*0.
   # Fake temperature matrix as there is no temperature on the Gill
   t_degC = np.zeros(len(u))*0.
   
   for j in range(len(u)):
      if u[j] != '':
         u_ms[j] = float(u[j])
         v_ms[j] = float(v[j])
         w_ms[j] = float(w[j])
      else:
         u_ms[j] = float('NaN')
         v_ms[j] = float('NaN')
         w_ms[j] = float('NaN')
         t_degC[j] = float('NaN')

   return u_ms, v_ms, w_ms, t_degC, datetimeList
            


