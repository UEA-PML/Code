#!/usr/bin/env python
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
import pylab 
import math
import argparse

rc('text',usetex=True)

# Function to create a moving average
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Function to allow three (or more) axes on a plot
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)


def toBinary(n, minb, maxb):
    return ''.join(str(1 & int(n) >> i) for i in range(64)[::-1])

def convert_as_needed(ts):
    dat_time = ''
    try:
        # parse strings
        dat_time = datetime.strptime(' '.join(ts), '%Y-%m-%d %H:%M:%S.%f')
    except:
        # assume doesn't have decimal fraction of a second
        print "error with timestamp string"
        ts = ts + ".5"
        dat_time = datetime.strptime(' '.join(ts), '%Y-%m-%d %H:%M:%S.%f')

    return dat_time


# Function to plot the data in the CR6
# LPMS, sonic, Licor
#"TIMESTAMP","SonicX","SonicY","SonicZ","SonicT","RotX","RotY","RotZ","AccX","AccY","AccZ","LicorCO2D","LicorH2OD","LicorTemp","LicorPress","LicorDiag"
def plot_CR6(outdir,ancillary):
   deg2rad = math.pi/180.
   rad2deg = 180./math.pi

   with open("/users/autoflux/tmp/latest_autoflux_CR6.txt") as f:
      CR6_data = f.readlines()

   tt = [row.split(',')[0:1] for row in CR6_data]
   x = [row.split(',')[1] for row in CR6_data]
   y = [row.split(',')[2] for row in CR6_data]
   z = [row.split(',')[3] for row in CR6_data]
   t = [row.split(',')[4] for row in CR6_data]
   rotx = [row.split(',')[5] for row in CR6_data]
   roty = [row.split(',')[6] for row in CR6_data]
   rotz = [row.split(',')[7] for row in CR6_data]
   accelx = [row.split(',')[8] for row in CR6_data] 
   accely = [row.split(',')[9] for row in CR6_data] 
   accelz = [row.split(',')[10] for row in CR6_data]
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

   lasttimestamp = datetime.strptime(''.join(tt[j]),'%Y-%m-%d %H:%M:%S.%f')
   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   fig.subplots_adjust(right=0.75)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,fco2_ppm,'r', label='CO$_{2}$')
   ax1.legend(loc=2)
   ax1.set_ylabel("CO$_{2}$ (ppm)")
   pylab.ylim([300,420])

   ax2 = ax1.twinx()
   ax2.xaxis.set_major_formatter(hfmt)
   ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax2.plot(datetimeList,fdiag,'b', label='Sig')
   ax2.set_ylabel(r"Signal Strength (\%)")
   ax2.legend(loc=1)
   pylab.ylim([0,100])
   
   ax3 = ax1.twinx()
   ax3.spines["right"].set_position(("axes", 1.2))
   make_patch_spines_invisible(ax3)   
   ax3.spines["right"].set_visible(True)
   
   ax3.xaxis.set_major_formatter(hfmt)
   ax3.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax3.plot(datetimeList,fh2o_ppm,'g', label='H$_{2}$O')
   ax3.set_ylabel("H$_{2}$O (ppm)")
   ax3.legend(loc=[0.815,0.82])
   pylab.ylim([30.,15000.])
   
   pylab.savefig(outdir+outfile_stem+'licor.png', dpi=100)

   if ancillary:
      print "Plotting Licor ancillary data"
      mave = 60
      fig = plt.figure()
      fig.set_size_inches(7.0,4.5)

      ax1 = fig.add_subplot(111)
   
      ax1.set_title("JCR ECCO2Flux - Preliminary Data")
      
      ax1.set_ylabel("Atmospheric pressure (atm)")
      pylab.ylim([960.,1050.])

      hfmt = dates.DateFormatter('%H:%M \n %d/%m')
      ax1.xaxis.set_major_formatter(hfmt)
      ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      plt.tick_params(axis='both', which='major', labelsize=10)   

      ax1.plot(datetimeList,fpressure,'b', label='Atmospheric pressure (mb)')
      ax1.legend(loc=2)

      ax2 = ax1.twinx()
      
      ax2.xaxis.set_major_formatter(hfmt)
      ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      ax2.plot(datetimeList,ftemp,'r--', label='Air temperature (degC)')
 
      ax2.set_ylabel("Temperature ($^\circ$C)")
      pylab.ylim([-10,30])
      ax2.legend(loc=1)
   
      pylab.savefig(outdir+outfile_stem+'licor_ancillary.png', dpi=100)

   ## LPMS ##
   frotx_rads = np.zeros(len(rotx))*0.
   froty_rads = np.zeros(len(roty))*0.
   frotz_rads = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_rads[j] = rad2deg*float(rotx[j])/1000.
      froty_rads[j] = rad2deg*float(roty[j])/1000.
      frotz_rads[j] = rad2deg*float(rotz[j])/1000.
      faccelx_g[j] = float(accelx[j])/1000.
      faccely_g[j] = float(accely[j])/1000.
      faccelz_g[j] = float(accelz[j])/1000.      

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Angular speed (deg s$^{-1}$)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,frotx_rads,'b', label='x')
   ax1.plot(datetimeList,froty_rads,'r', label='y')
   ax1.plot(datetimeList,frotz_rads,'g', label='z')
   
   pylab.ylim([-2.0,2.0])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'rot_motion.png', dpi=100)

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Acceleration (g)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,faccelx_g,'b', label='x')
   ax1.plot(datetimeList,faccely_g,'r', label='y')
   ax1.plot(datetimeList,faccelz_g,'g', label='z')
   
   pylab.ylim([-1.3,0.1])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'accel_motion.png', dpi=100)

# Function to plot the CO2 and CH4 data from the Picarro
def plot_co2_ch4(outdir, ancillary, picarro_hdd):

   if picarro_hdd:
      print "Using alternative Picarro HDD source"
      with open("/users/autoflux/tmp/latest_autoflux_picarro_picarro.txt") as f:
         picarro_data = f.readlines()
      
      tt = [row.split(' ')[0:2] for row in picarro_data]
      pump_pressure_torr = [row.split(' ')[8] for row in picarro_data]
      co2_temperature = [row.split(' ')[9] for row in picarro_data]
      switch_state = [row.split(' ')[16] for row in picarro_data]

      ch4 = [row.split(' ')[18] for row in picarro_data]
      co2_air = [row.split(' ')[19] for row in picarro_data]
      co2_dry = [row.split(' ')[20] for row in picarro_data]
      h2o = [row.split(' ')[21] for row in picarro_data]

      datetimeList = [datetime.strptime(' '.join(row),'%Y-%m-%d %H:%M:%S.%f') for row in tt]
      
   else:
      with open("/users/autoflux/tmp/latest_autoflux_picarro.txt") as f:
         picarro_data = f.readlines()
   
      tt = [row.split(' ')[0:5] for row in picarro_data]
      seconds = [row.split(' ')[5] for row in picarro_data]
      pump_pressure_torr = [row.split(' ')[6] for row in picarro_data]
      co2_temperature = [row.split(' ')[7] for row in picarro_data]
      switch_state = [row.split(' ')[8] for row in picarro_data]
   
      ch4 = [row.split(' ')[9] for row in picarro_data]
      co2_air = [row.split(' ')[10] for row in picarro_data]
      co2_dry = [row.split(' ')[11] for row in picarro_data]
      h2o = [row.split(' ')[12] for row in picarro_data]

      datetimeList = [datetime.strptime(' '.join(row),'%a %b %d %H:%M:%S.%f %Y') for row in tt]
   
   
   fpump_pressure_atm_co2 = np.zeros(len(pump_pressure_torr))*0.
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
      fpump_pressure_atm_co2[j] = float(pump_pressure_torr[j])/760.
      fco2_temperature[j] = float(co2_temperature[j])
      fswitch_state[j] = float(switch_state[j]) - 32.0

   if picarro_hdd:
      lasttimestamp = datetime.strptime(' '.join(tt[j]),'%Y-%m-%d %H:%M:%S.%f')
   else:
      lasttimestamp = datetime.strptime(' '.join(tt[j]),'%a %b %d %H:%M:%S.%f %Y')

   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'
   
   # CH4 and CO2
   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)

   fig.subplots_adjust(right=0.75)

   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   #ax1.plot(datetimeList,fco2_air,'r', label='CO$_{2}$ air')
   ax1.plot(datetimeList,fco2_dry,'r--', label='CO$_{2}$ dry')
   ax1.legend(loc=2)
   ax1.set_ylabel("CO$_{2}$ (ppm)")
   pylab.ylim([395,410])

   ax2 = ax1.twinx()
   ax2.xaxis.set_major_formatter(hfmt)
   ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax2.plot(datetimeList,fch4,'b', label='CH$_{4}$')
   ax2.set_ylabel("CH$_{4}$ (ppm)")
   ax2.legend(loc=1)
   pylab.ylim([1.6,2.2])
   
   ax3 = ax1.twinx()
   ax3.spines["right"].set_position(("axes", 1.2))
   make_patch_spines_invisible(ax3)   
   ax3.spines["right"].set_visible(True)
   
   ax3.xaxis.set_major_formatter(hfmt)
   ax3.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax3.plot(datetimeList,fh2o,'g', label='H$_{2}$O')
   ax3.plot(datetimeList,fswitch_state,'k', label='switch')
   ax3.set_ylabel(r"H$_{2}$O (\%)")
   ax3.legend(loc=[0.80,0.76])
   pylab.ylim([0.,0.5])
   
   pylab.savefig(outdir+outfile_stem+'co2_ch4.png', dpi=100)

   if ancillary:
      print "Plotting ancillary data"
      mave = 60
      fig = plt.figure()
      fig.set_size_inches(7.0,4.5)

      ax1 = fig.add_subplot(111)
   
      ax1.set_title("JCR ECCO2Flux - Preliminary Data")
      
      ax1.set_ylabel("pump pressure (atm)")
      pylab.ylim([0.,1.1])

      hfmt = dates.DateFormatter('%H:%M \n %d/%m')
      ax1.xaxis.set_major_formatter(hfmt)
      ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      plt.tick_params(axis='both', which='major', labelsize=10)   

      ax1.plot(datetimeList[mave-1:len(datetimeList)],moving_average(fpump_pressure_atm_co2,n=mave),'b', label='CO$_{2}$ pump press.')
      ax1.legend(loc=2)

      ax2 = ax1.twinx()
      
      ax2.xaxis.set_major_formatter(hfmt)
      ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      ax2.plot(datetimeList[mave-1:len(datetimeList)],moving_average(fco2_temperature,n=mave),'b--', label='CO$_{2}$ temperature')
 
      ax2.set_ylabel("Temperature ($^\circ$C)")
      pylab.ylim([0,70])
      ax2.legend(loc=1)
   
      pylab.savefig(outdir+outfile_stem+'pump_pressure.png', dpi=100)

# Function to plot the Licor data
def plot_licor(outdir, ancillary):
   with open("/users/autoflux/tmp/latest_autoflux_licor.txt") as f:
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
      fpressure[j] = float(pressure[j])*10. #/101.325
      fco2_ppm[j] = 1e+3*fco2_density[j]/(fPRT[j] - fh2o_density[j]/1000.)
      fh2o_ppm[j] = 1e+3*fh2o_density[j]/fPRT[j]
   
   # Filename created based on the latest date available
   lasttimestamp = datetime.strptime(' '.join(tt[j]),'%a %b %d %H:%M:%S.%f %Y')
   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   fig.subplots_adjust(right=0.75)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,fco2_ppm,'r', label='CO$_{2}$')
   ax1.legend(loc=2)
   ax1.set_ylabel("CO$_{2}$ (ppm)")
   pylab.ylim([300,420])

   ax2 = ax1.twinx()
   ax2.xaxis.set_major_formatter(hfmt)
   ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax2.plot(datetimeList,fdiag,'b', label='Sig')
   ax2.set_ylabel(r"Signal Strength (\%)")
   ax2.legend(loc=1)
   pylab.ylim([0,100])
   
   ax3 = ax1.twinx()
   ax3.spines["right"].set_position(("axes", 1.2))
   make_patch_spines_invisible(ax3)   
   ax3.spines["right"].set_visible(True)
   
   ax3.xaxis.set_major_formatter(hfmt)
   ax3.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,3)))
   ax3.plot(datetimeList,fh2o_ppm,'g', label='H$_{2}$O')
   ax3.set_ylabel("H$_{2}$O (ppm)")
   ax3.legend(loc=[0.815,0.82])
   pylab.ylim([30.,15000.])
   
   pylab.savefig(outdir+outfile_stem+'licor.png', dpi=100)
   
   if ancillary:
      print "Plotting Licor ancillary data"
      mave = 60
      fig = plt.figure()
      fig.set_size_inches(7.0,4.5)

      ax1 = fig.add_subplot(111)
   
      ax1.set_title("JCR ECCO2Flux - Preliminary Data")
      
      ax1.set_ylabel("Atmospheric pressure (atm)")
      pylab.ylim([960.,1050.])

      hfmt = dates.DateFormatter('%H:%M \n %d/%m')
      ax1.xaxis.set_major_formatter(hfmt)
      ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      plt.tick_params(axis='both', which='major', labelsize=10)   

      ax1.plot(datetimeList,fpressure,'b', label='Atmospheric pressure (mb)')
      ax1.legend(loc=2)

      ax2 = ax1.twinx()
      
      ax2.xaxis.set_major_formatter(hfmt)
      ax2.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
      ax2.plot(datetimeList,ftemp,'r--', label='Air temperature (degC)')
 
      ax2.set_ylabel("Temperature ($^\circ$C)")
      pylab.ylim([-10,30])
      ax2.legend(loc=1)

      pylab.savefig(outdir+outfile_stem+'licor_ancillary.png', dpi=100)

def plot_motion(outdir):
   with open("/users/autoflux/tmp/latest_autoflux_LPMS.txt") as f:
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

   frotx_rads = np.zeros(len(rotx))*0.
   froty_rads = np.zeros(len(roty))*0.
   frotz_rads = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_rads[j] = rad2deg*float(rotx[j])/1000.
      froty_rads[j] = rad2deg*float(roty[j])/1000.
      frotz_rads[j] = rad2deg*float(rotz[j])/1000.
      faccelx_g[j] = float(accelx[j])/1000.
      faccely_g[j] = float(accely[j])/1000.
      faccelz_g[j] = float(accelz[j])/1000.      

   lasttimestamp = datetime.strptime(' '.join(tt[j]),'%Y %m %d %H %M %S.%f')
   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Angular speed (deg s$^{-1}$)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,frotx_rads,'b', label='x')
   ax1.plot(datetimeList,froty_rads,'r', label='y')
   ax1.plot(datetimeList,frotz_rads,'g', label='z')
   
   pylab.ylim([-2.0,2.0])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'rot_motion.png', dpi=100)

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Acceleration (g)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,faccelx_g,'b', label='x')
   ax1.plot(datetimeList,faccely_g,'r', label='y')
   ax1.plot(datetimeList,faccelz_g,'g', label='z')
   
   pylab.ylim([-1.3,0.1])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'accel_motion.png', dpi=100)

def plot_sysdon(outdir):
   with open("/users/autoflux/tmp/latest_autoflux_SysDon.txt") as f:
      motion_data = f.readlines()
   
   deg2rad = math.pi/180.
   
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

   frotx_rads = np.zeros(len(rotx))*0.
   froty_rads = np.zeros(len(roty))*0.
   frotz_rads = np.zeros(len(rotz))*0.
   faccelx_g = np.zeros(len(accelx))*0.
   faccely_g = np.zeros(len(accely))*0.
   faccelz_g = np.zeros(len(accelz))*0.
   
   for j in range(len(rotx)):
      frotx_rads[j] = -1.0*(((float(rotx[j])*5./4095.)-2.5)/0.027)
      froty_rads[j] = -1.0*(((float(roty[j])*5./4095.)-2.5)/0.027)
      frotz_rads[j] = -1.0*(((float(rotz[j])*5./4095.)-2.5)/0.027)
      faccelx_g[j] = -1.0*(((float(accelx[j])*5./4095.)-2.5)/0.75)
      faccely_g[j] = -1.0*(((float(accely[j])*5./4095.)-2.5)/0.75)
      faccelz_g[j] = -1.0*(((float(accelz[j])*5./4095.)-2.5)/0.75)      

   lasttimestamp = datetime.strptime(''.join(tt[j]),'%Y-%m-%d %H:%M:%S.%f')
   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Angular speed (deg s$^{-1}$)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,frotx_rads,'b', label='x')
   ax1.plot(datetimeList,froty_rads,'r', label='y')
   ax1.plot(datetimeList,frotz_rads,'g', label='z')
   
   pylab.ylim([-2.0,2.0])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'SysDon_rot_motion.png', dpi=100)

   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Acceleration (g)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,faccelx_g,'b', label='x')
   ax1.plot(datetimeList,faccely_g,'r', label='y')
   ax1.plot(datetimeList,faccelz_g,'g', label='z')
   
   pylab.ylim([-1.3,0.1])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'SysDon_accel_motion.png', dpi=100)


def plot_wind(outdir):
   with open("/users/autoflux/tmp/latest_autoflux_wind.txt") as f:
      wind_data = f.readlines()
   
   tt = [row.split(' ')[0:5] for row in wind_data]
   x = [row.split(' ')[5] for row in wind_data]
   y = [row.split(' ')[6] for row in wind_data]
   z = [row.split(' ')[7] for row in wind_data]
   t = [row.split(' ')[8] for row in wind_data]
   
   datetimeList = [datetime.strptime(' '.join(row),'%a %b %d %H:%M:%S.%f %Y') for row in tt]
   
   u_ms = np.zeros(len(x))*0.
   v_ms = np.zeros(len(y))*0.
   w_ms = np.zeros(len(z))*0.
   
   wind_spd = np.zeros(len(x))*0.
   wind_dir = np.zeros(len(x))*0.
   
   for j in range(len(x)):
      if x[j] != '':
         u_ms[j] = float(x[j])/100.
         v_ms[j] = float(y[j])/100.
         w_ms[j] = float(z[j])/100.
         wind_spd[j] = math.sqrt(u_ms[j]**2 + v_ms[j]**2)
         wind_dir[j] = math.atan2(v_ms[j],-u_ms[j])*(180./math.pi)
         if wind_dir[j] < 0.:
            wind_dir[j] += 360.
      else:
         wind_spd[j] = float('NaN')
         wind_dir[j] = float('NaN')

   #datetime.strptime(' '.join(tt[j]),'%a %b %d %H:%M:%S.%f %Y'),strFormat="%Y-%m-%d")    
   #timeStamp = time.mktime(time.strptime(str_date,strFormat))
   lasttimestamp = datetime.strptime(' '.join(tt[j]),'%a %b %d %H:%M:%S.%f %Y')
   isodate = lasttimestamp.isoformat()
   outfile_stem = isodate[0:10]+'_'

   # Wind speed and direction
   fig = plt.figure()
   fig.set_size_inches(7.0,4.5)
   ax1 = fig.add_subplot(111)
   ax1.set_title("JCR ECCO2Flux - Preliminary Data")
   
   hfmt = dates.DateFormatter('%H:%M \n %d/%m')

   ax1.set_ylabel("Wind Velocity (ms$^{-1}$)")
   ax1.xaxis.set_major_formatter(hfmt)
   ax1.xaxis.set_major_locator(dates.HourLocator(byhour=range(0,24,2)))
   plt.tick_params(axis='both', which='major', labelsize=10)   
   ax1.plot(datetimeList,u_ms,'b', label='u')
   ax1.plot(datetimeList,v_ms,'r', label='v')
   ax1.plot(datetimeList,w_ms,'g', label='w')
   
   pylab.ylim([-20,20])
   ax1.legend(loc=2)

   pylab.savefig(outdir+outfile_stem+'wind.png', dpi=100)
   #pylab.savefig(outdir+'wind.png', dpi=300)
      
def main():

   parser = argparse.ArgumentParser(
      description=__doc__,
      formatter_class=argparse.RawDescriptionHelpFormatter
   )

   parser.add_argument("-w", "--wind", action='store_true', default=False, help="Plot the wind data.")
   parser.add_argument("-o", "--outdir", type=str, default='', help="Output directory for png imagery")
   parser.add_argument("-p", "--picarro", action='store_true', default=False, help="Plot the Picarro data.")
   parser.add_argument("-pd", "--picarro_hdd", action='store_true', default=False, help="Plot the Picarro data from alternative Picarro HDD source.")
   parser.add_argument("-a", "--ancillary", action='store_true', default=False, help="Plot the ancillary data.")
   parser.add_argument("-l", "--licor", action='store_true', default=False, help="Plot the Licor data.")
   parser.add_argument("-m", "--motion", action='store_true', default=False, help="Plot the LPMS motion data.")
   parser.add_argument("-s", "--sysdon", action='store_true', default=False, help="Plot the SysDon motion data.")
   parser.add_argument("-c", "--CR6", action='store_true', default=False, help="Plot the CR6 data streams.")
   args = parser.parse_args()
   
   if args.wind:
      print "Producing wind plot"
      plot_wind(args.outdir)

   if args.picarro:
      print "Producing Picarro CO2 and CH4 plot"
      plot_co2_ch4(args.outdir, args.ancillary, args.picarro_hdd)

   if args.picarro_hdd:
      print "Producing Picarro CO2 and CH4 plot using alternative Picarro data source"
      plot_co2_ch4(args.outdir, args.ancillary, args.picarro_hdd)

   if args.licor:
      print "Producing Licor plot"
      plot_licor(args.outdir, args.ancillary)
   
   if args.motion:
      print "Producing motion plot from LPMS"
      plot_motion(args.outdir)

   if args.sysdon:
      print "Producing motion plot from SysDon"
      plot_sysdon(args.outdir)
   
   if args.CR6:
      print "Producing CR6 plots"
      plot_CR6(args.outdir, args.ancillary)
         
   sys.exit()
      
if __name__=='__main__':
   main()
   
