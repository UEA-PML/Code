#export PATH="/Users/MXY/Python2/anaconda2/bin:$PATH"
TEST SOME CHANGE
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
from numpy import copy, array, nanmean
from numpy import cos, sin
import pylab 
import math
import argparse
import pandas as pd
from scipy import signal
import argparse
import csv
import glob
import os

# set current directory for the eddy covariance python scripts
os.chdir('/Users/mingxi/Dropbox/ecco2flux')
#os.chdir('/Users/MXY/Documents/Students/Yuanxu Dong/Python ECco2flux scripts JCR')
from EC_functions import *
#from read_instruments import *
from RunMotionCorrection import *

# set current directory for the raw datafiles
os.chdir('/Users/mingxi/Desktop/WinddataJCRApr2019/')
#os.chdir('/Users/MXY/Desktop/WinddataJCRApr2019/')

def RunMotCorrection():
	######## Definition of fixed parameters ##########
	deg2rad = math.pi/180.    # conversion of degrees to radians
	rad2deg = 180./math.pi    # conversion of radians to degrees
	g  = 9.80665              # gravitational acceleration (m/s2) - conventional standard value
	knots2ms = 0.51444444     # conversion of knots to m/s

	# Define position vector for distance between MotionPak and sonic volume.
	# Units are metres, x positive forward, y positive to port, z positive up.
	# These values are for the system on the JCR
	Rx = 0.17   #Metek Sonic is 17 cm in front of MotionPak
	Ry = 0.  
	Rz = 0.725  #Metek Sonic is 72.5 cm above MotionPak
	alpha = -0.3 #Metek Sonic rotated 0.3 degrees towards port side relative to the motion sensor
	beta = 1.2  #Metek Sonic rotated 1.2 degrees towards the stern relative to the motion sensor

	fsam = 10   # All data to be interpolated onto a 10 Hz grid
	# names of the data files to be loaded, sorted by number
	filenames = sorted(glob.glob('WindMotion*.txt'))

	# make 1-pt array for summary data. 
	# this has because concatenation function vstack doesn't like to append to empty array
	FileTime20_save = np.zeros(1)*0. 
 	U20_save = np.zeros(1)*0. 
 	Uvar20_save = np.zeros(1)*0. 
 	V20_save = np.zeros(1)*0. 
 	Vvar20_save = np.zeros(1)*0. 
 	W20_save = np.zeros(1)*0. 
 	Wvar20_save = np.zeros(1)*0. 
 	Uplat20_save = np.zeros(1)*0. 
 	UplatVar20_save = np.zeros(1)*0. 
 	Vplat20_save = np.zeros(1)*0. 
 	VplatVar20_save = np.zeros(1)*0. 
 	Wplat20_save = np.zeros(1)*0. 
 	WplatVar20_save = np.zeros(1)*0. 
 	U20true_save = np.zeros(1)*0. 
 	UtrueVar20_save = np.zeros(1)*0. 
 	RelWdir20_save = np.zeros(1)*0. 
 	Ustar20_cov_save = np.zeros(1)*0. 
 	TrueWdir20_save = np.zeros(1)*0. 
	SOS20_save = np.zeros(1)*0. 
 	Tilt20_save = np.zeros(1)*0. 	

	for ii in range(len(filenames)):
		file2load = filenames[ii]
		with open(file2load, "r") as f:

			#contents=f.readlines()
			# skip headerline
			contents=f.readlines()[1:]

		# data initially read in as strings
		datesecs = [row.split('\t')[0] for row in contents]
		u = [row.split('\t')[1] for row in contents]
		v = [row.split('\t')[2] for row in contents]
		w = [row.split('\t')[3] for row in contents]
		sonicT = [row.split('\t')[4] for row in contents]
		x_accel = [row.split('\t')[5] for row in contents]
		y_accel = [row.split('\t')[6] for row in contents]
		z_accel = [row.split('\t')[7] for row in contents]
		x_rot = [row.split('\t')[8] for row in contents]
		y_rot = [row.split('\t')[9] for row in contents]
		z_rot = [row.split('\t')[10] for row in contents]
		smg_knots = [row.split('\t')[11] for row in contents]
		heading_deg = [row.split('\t')[12] for row in contents]

		# create floating numerical equivalent placeholders and set to 0
		fdatesecs = np.zeros(len(datesecs))*0.
		fu = np.zeros(len(u))*0.
		fv = np.zeros(len(v))*0.
		fw = np.zeros(len(w))*0.
		fsonicT = np.zeros(len(sonicT))*0.
		fx_accel = np.zeros(len(x_accel))*0.
		fy_accel = np.zeros(len(y_accel))*0.
		fz_accel = np.zeros(len(z_accel))*0.
		fx_rot = np.zeros(len(x_rot))*0.
		fy_rot = np.zeros(len(y_rot))*0.
		fz_rot = np.zeros(len(z_rot))*0.
		fsmg_knots = np.zeros(len(smg_knots))*0.
		fheading_deg = np.zeros(len(heading_deg))*0.

		# fill in the values   
		for j in range(len(datesecs)):
			fdatesecs[j] = float(datesecs[j])
			fu[j] = float(u[j])
			fv[j] = float(v[j])
			fw[j] = float(w[j])
			fsonicT[j] = float(sonicT[j])
			fx_accel[j] = float(x_accel[j])
			fy_accel[j] = float(y_accel[j]) 
			fz_accel[j] = float(z_accel[j]) 
			fx_rot[j] = float(x_rot[j]) 
			fy_rot[j] = float(y_rot[j]) 
			fz_rot[j] = float(z_rot[j]) 
			fsmg_knots[j] = float(smg_knots[j]) 
			fheading_deg[j] = float(heading_deg[j]) 

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
		ratex = -fy_rot * deg2rad
		ratey =  fx_rot * deg2rad
		ratez =  fz_rot * deg2rad
		accx  = -fy_accel * g 
		accy  =  fx_accel * g
		accz  =  fz_accel * g

		# Convert Speed Of Ship from knots to m/s
		SOS = fsmg_knots * knots2ms   

		# Get ship's heading in radians for use in yaw angle and true wind calculation.
		heading = np.unwrap(fheading_deg * deg2rad) 

		# correct Metek data for angular offset of installation
		polarity = 1	# Metek
		# correct Metek data for angular offset of installation
		u_ms, v_ms, w_ms = angular_offset(fu, fv, fw, alpha*deg2rad, beta*deg2rad, polarity)

		# set up arrays for motion correction
		#son: Nx4 array of time and sonic wind velocities.
		#    N is typically 36000 for 10 Hz data.  
		son = np.zeros((len(fu),4))*0.
		for j in range(len(fu)):
		    son[j,0] = fdatesecs[j]
		    son[j,1] = u_ms[j]
		    son[j,2] = v_ms[j]
		    son[j,3] = w_ms[j]

		#mot: Nx7 array of time, linear acceleration and angular rate variables.
		mot = np.zeros((len(fu),7))*0.
		for j in range(len(fu)):
			mot[j,0] = fdatesecs[j]
			mot[j,1] = accx[j]
			mot[j,2] = accy[j]
			mot[j,3] = accz[j]
			mot[j,4] = ratex[j]
			mot[j,5] = ratey[j]
			mot[j,6] = ratez[j]

		# displacements
		sens_disp = np.zeros(3)*0.
		sens_disp[0] = Rx
		sens_disp[1] = Ry
		sens_disp[2] = Rz
	
		# perform motion correction
		# uvw is the motion-corrected winds following Edson et al. 1998
		acc, uvw, uvwplat, xyzplat = motcorr(son,mot,heading,sens_disp)

		uvw_mc = uvw
		#now reverse the sign of u to bring into meteorological convention
		uvw_mc[:,0] = -uvw[:,0]
			
		#### additional Decorrelation correction between winds and motion ####
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
	   	########

		# Correct winds for streamline distortion by ship superstructure.
		# Streamline corrected uvw (NOT corrected for ship speed) will be used
		# for flux calculations.  
		uvw_str, tilt, azm, Umean, Vmean, Wmean = dbl_rot(uvw_mc)

		# Compute true wind speed, relative and true wind directions
		#In vectors: Relative Wind = True Wind + ship velocity over water.
		#    Convert to rectangular coordinates, then 
		#    RelWind(X2,Y2) - SOS(X1,Y1) = TrueWind(X3,Y3)
		#    and then convert back to polar coordinates.
		# Utrue in this calculation is wind speed relative to the ocean
		#    surface. (i.e. accounts for ocean current)
		x1 = np.zeros(len(fu))*0.
		x2 = np.zeros(len(fu))*0.
		y1 = np.zeros(len(fu))*0.
		y2 = np.zeros(len(fu))*0.
		x3 = np.zeros(len(fu))*0.
		y3 = np.zeros(len(fu))*0.
		TrueWdir = np.zeros(len(fu))*0.
		Utrue = np.zeros(len(fu))*0.

		RelWdirGeo = np.zeros(len(fu))*0.
		RelWdir = np.zeros(len(fu))*0.
   
		for k in range(len(fu)):
			RelWdirGeo[k] = rad2deg*np.arctan2(uvw_mc[k, 1],uvw_mc[k, 0])
			if (RelWdirGeo[k] > 360.):
				RelWdirGeo[k] = RelWdirGeo[k] - 360.
			if (RelWdirGeo[k] < -360.):
				RelWdirGeo[k] = RelWdirGeo[k] + 360.
			RelWdir[k] = RelWdirGeo[k] - (fheading_deg[k]-360.)
			if (RelWdir[k] > 180.):
				RelWdir[k] = RelWdir[k] - 360. # +/- 180 deg range 
			if (RelWdir[k] < -180.):
				RelWdir[k] = RelWdir[k] + 360. # +/- 180 deg range 
       
			x1[k], y1[k] = pol2cart(fheading_deg[k],SOS[k])
			# uvw_str is streamline corrected uvw
			x2[k], y2[k] = pol2cart(RelWdirGeo[k],uvw_str[k,0])
			x3[k] = x2[k] - x1[k]
			y3[k] = y2[k] - y1[k]
			TrueWdir[k], Utrue[k] = cart2pol(x3[k],y3[k]) # cart2pol returns direction in degrees
			if (TrueWdir[k] > 360.): # check for overly positive angles
				TrueWdir[k] = TrueWdir[k] - 360.
			if (TrueWdir[k] < 0.): # check for negative angles
				TrueWdir[k] = TrueWdir[k] + 360.

		# save 10-Hz motion corrected data
		Data2save = np.zeros((len(fu),17))*0.
		for j in range(len(fu)):
			Data2save[j,0] = round(fdatesecs[j], 2)
			Data2save[j,1] = round(uvw_str[j][0], 4)
			Data2save[j,2] = round(uvw_str[j][1], 4)
			Data2save[j,3] = round(uvw_str[j][2], 4)
			Data2save[j,4] = round(RelWdir[j], 2)
			Data2save[j,5] = round(Utrue[j], 4)
			Data2save[j,6] = round(fheading_deg[j], 2)
			Data2save[j,7] = round(fsonicT[j], 3)
			Data2save[j,8] = round(acc[j][0], 4)
			Data2save[j,9] = round(acc[j][1], 4)
			Data2save[j,10] = round(acc[j][2], 4)
			Data2save[j,11] = round(uvwplat[j][0], 4)
			Data2save[j,12] = round(uvwplat[j][1], 4)
			Data2save[j,13] = round(uvwplat[j][2], 4)
			Data2save[j,14] = round(xyzplat[j][0], 4)
			Data2save[j,15] = round(xyzplat[j][1], 4)
			Data2save[j,16] = round(xyzplat[j][2], 4)
		
		savename = file2load[10:20] + 'motcorr.txt'

		with open(savename, 'w') as a:
			csv.writer(a, delimiter='	').writerows(Data2save)

		#########################################################
		# Calculate 20 min wind & motion statistics as well as u*

		# unwrap wind direction for calculations below
		TrueWdir = np.unwrap(TrueWdir*deg2rad)*rad2deg
		
		# define starting and ending ponts for 20 min segment
		startPt = 1;
		endPt = fsam*60*20;  # number of points in 20 min segment, here = 1st end pt
		lastPt = len(fu);	# very last point of hour
				
		if (endPt > lastPt):
			endPt = lastPt;

		# number of 20-min segments within an hour
		# for partial hour, lastz<3
		lastz = round(float(lastPt)/(fsam*60*20)); 
		if (lastz > 3):
			lastz = 3;
			  		
		lastz = int(lastz)  
		FileTime20 = np.zeros(lastz)*0.  
		U20 = np.zeros(lastz)*0. 
		Uvar20 = np.zeros(lastz)*0. 
		V20 = np.zeros(lastz)*0. 
		Vvar20 = np.zeros(lastz)*0. 
		W20 = np.zeros(lastz)*0. 
		Wvar20 = np.zeros(lastz)*0. 
		Uplat20 = np.zeros(lastz)*0. 
		UplatVar20 = np.zeros(lastz)*0. 
		Vplat20 = np.zeros(lastz)*0. 
		VplatVar20 = np.zeros(lastz)*0. 
		Wplat20 = np.zeros(lastz)*0. 
		WplatVar20 = np.zeros(lastz)*0. 
		U20true = np.zeros(lastz)*0. 
		UtrueVar20 = np.zeros(lastz)*0. 
		RelWdir20 = np.zeros(lastz)*0. 
		Ustar20_cov = np.zeros(lastz)*0. 
		TrueWdir20 = np.zeros(lastz)*0. 
		SOS20 = np.zeros(lastz)*0. 
		Tilt20 = np.zeros(lastz)*0.
				
		# compute stats
		for z in range(lastz):
			FileTime20[z] = np.mean(fdatesecs[startPt:endPt+1]);
			U20[z] = np.mean(uvw_str[startPt:endPt+1, 0]);
			Uvar20[z] = np.var(uvw_str[startPt:endPt+1, 0]);
			V20[z] = np.mean(uvw_str[startPt:endPt+1, 1]);
			Vvar20[z] = np.var(uvw_str[startPt:endPt+1, 1]);
			W20[z] = np.mean(uvw_str[startPt:endPt+1, 2]);
			Wvar20[z] = np.var(uvw_str[startPt:endPt+1, 2]);
			Uplat20[z] = np.mean(uvwplat[startPt:endPt+1, 0]);
			UplatVar20[z] = np.var(uvwplat[startPt:endPt+1, 0]);
			Vplat20[z] = np.mean(uvwplat[startPt:endPt+1, 1]);
			VplatVar20[z] = np.var(uvwplat[startPt:endPt+1, 1]);
			Wplat20[z] = np.mean(uvwplat[startPt:endPt+1, 2]);
			WplatVar20[z] = np.var(uvwplat[startPt:endPt+1, 2]);
			U20true[z] = np.mean(Utrue[startPt:endPt+1]);
			UtrueVar20[z] = np.var(Utrue[startPt:endPt+1]);
			RelWdir20[z] = np.mean(RelWdir[startPt:endPt+1]);
			TrueWdir20[z] = np.mean(TrueWdir[startPt:endPt+1]);
			TrueWdir20[z] = np.mod(TrueWdir20[z], 360.)
			if (RelWdir20[z] > 180.):
				RelWdir20[z] = RelWdir20[z]-360.
			if (RelWdir20[z] < -180.):
				RelWdir20[z] = RelWdir20[z]+360.
			if (TrueWdir20[z] > 360.):
				TrueWdir20[z] = TrueWdir20[z]-360.
			if (TrueWdir20[z] < -360.):
				TrueWdir20[z] = TrueWdir20[z]+360.
			Ustar20_cov[z] = np.sqrt(np.absolute(np.mean((uvw_str[startPt:endPt+1, 0]-U20[z])*(uvw_str[startPt:endPt+1, 2]-W20[z]))));
			SOS20[z] = np.mean(SOS[startPt:endPt+1]);
			Tilt20[z] = tilt;
			
			startPt = startPt + (fsam*60*20);
			endPt = startPt + (fsam*60*20)-1;
			if (endPt > lastPt):
				endPt = lastPt;

 		# concatenate results
 		FileTime20_save = np.hstack((FileTime20_save,FileTime20))
 		U20_save = np.hstack((U20_save,U20))
 		Uvar20_save = np.hstack((Uvar20_save,Uvar20))
 		V20_save = np.hstack((V20_save,V20))
 		Vvar20_save = np.hstack((Vvar20_save,Vvar20))
 		W20_save = np.hstack((W20_save,W20))
 		Wvar20_save = np.hstack((Wvar20_save,Wvar20))
 		Uplat20_save = np.hstack((Uplat20_save,Uplat20))
 		UplatVar20_save = np.hstack((UplatVar20_save,UplatVar20))
 		Vplat20_save = np.hstack((Vplat20_save,Vplat20))
 		VplatVar20_save = np.hstack((VplatVar20_save,VplatVar20))
 		Wplat20_save = np.hstack((Wplat20_save,Wplat20))
 		WplatVar20_save = np.hstack((WplatVar20_save,WplatVar20))
 		U20true_save = np.hstack((U20true_save,U20true))
 		UtrueVar20_save = np.hstack((UtrueVar20_save,UtrueVar20))
 		RelWdir20_save = np.hstack((RelWdir20_save,RelWdir20))
 		Ustar20_cov_save = np.hstack((Ustar20_cov_save,Ustar20_cov))
 		TrueWdir20_save = np.hstack((TrueWdir20_save,TrueWdir20))
 		SOS20_save = np.hstack((SOS20_save,SOS20))
 		Tilt20_save = np.hstack((Tilt20_save,Tilt20))		

	# save 20-minute summary
	# length -1 here (and index j+1 below) to account for the fact that the first empty (0) is not saved
	SumData2save = np.zeros((len(FileTime20_save)-1,20))*0.

	for j in range(len(FileTime20_save)-1):
		SumData2save[j,0] = round(FileTime20_save[j+1], 0)
		SumData2save[j,1] = round(U20_save[j+1], 3)
		SumData2save[j,2] = round(Uvar20_save[j+1], 3)
		SumData2save[j,3] = round(V20_save[j+1], 3)
		SumData2save[j,4] = round(Vvar20_save[j+1], 3)
		SumData2save[j,5] = round(W20_save[j+1], 3)
		SumData2save[j,6] = round(Wvar20_save[j+1], 3)
		SumData2save[j,7] = round(Uplat20_save[j+1], 4)
		SumData2save[j,8] = round(UplatVar20_save[j+1], 4)
		SumData2save[j,9] = round(Vplat20_save[j+1], 4)
		SumData2save[j,10] = round(VplatVar20_save[j+1], 4)
		SumData2save[j,11] = round(Wplat20_save[j+1], 4)
		SumData2save[j,12] = round(WplatVar20_save[j+1], 4)
		SumData2save[j,13] = round(U20true_save[j+1], 4)
		SumData2save[j,14] = round(UtrueVar20_save[j+1], 3)
		SumData2save[j,15] = round(RelWdir20_save[j+1], 2)
		SumData2save[j,16] = round(Ustar20_cov_save[j+1], 4)
		SumData2save[j,17] = round(TrueWdir20_save[j+1], 2)
		SumData2save[j,18] = round(SOS20_save[j+1], 2)
		SumData2save[j,19] = round(Tilt20_save[j+1], 2)

	savename = 'SumStat20min.txt'

	with open(savename, 'w') as a:
		csv.writer(a, delimiter='	').writerows(SumData2save)

	return ;
RunMotCorrection()
