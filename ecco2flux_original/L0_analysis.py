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

   freq_data_matrix_gill = pd.read_pickle('output/L0/Gill/20181018-195958_L0.pkl')
      
   U_gill = copy(freq_data_matrix_gill.values[:,0])
   V_gill = copy(freq_data_matrix_gill.values[:,1])
   W_gill = copy(freq_data_matrix_gill.values[:,2])
   timestamps = freq_data_matrix_gill.axes
   timestamps = timestamps[0]

   U_gill_ts = pd.Series(U_gill, index=timestamps)
   U_gill_ts = U_gill_ts.resample('1s').mean()
   V_gill_ts = pd.Series(V_gill, index=timestamps)
   V_gill_ts = V_gill_ts.resample('1s').mean()
   W_gill_ts = pd.Series(W_gill, index=timestamps)
   W_gill_ts = W_gill_ts.resample('1s').mean()

   freq_data_matrix_metek = pd.read_pickle('output/L0/Metek/20181018-195958_L0.pkl')
      
   U_metek = copy(freq_data_matrix_metek.values[:,0])
   V_metek = copy(freq_data_matrix_metek.values[:,1])
   W_metek = copy(freq_data_matrix_metek.values[:,2])
   timestamps = freq_data_matrix_metek.axes
   timestamps = timestamps[0]

   U_metek_ts = pd.Series(U_metek, index=timestamps)
   U_metek_ts = U_metek_ts.resample('1s').mean()
   V_metek_ts = pd.Series(V_metek, index=timestamps)
   V_metek_ts = V_metek_ts.resample('1s').mean()
   W_metek_ts = pd.Series(W_metek, index=timestamps)
   W_metek_ts = W_metek_ts.resample('1s').mean()

   U_residuals = np.mean(U_metek_ts.values - U_gill_ts.values)
   V_residuals = np.mean(V_metek_ts.values - V_gill_ts.values)
   W_residuals = np.mean(W_metek_ts.values - W_gill_ts.values)
   
   print U_residuals
   print V_residuals
   print W_residuals
      
   fig = plt.figure()
   fig.set_size_inches(5,3)
   ax1 = fig.add_subplot(111)
   ax1.set_title("Metek \& Gill intercomparison")
   ax1.set_xlabel("Metek U(m/s)")
   ax1.set_ylabel("Gill U(m/s)")
   pylab.ylim([-10,10])
   pylab.xlim([-10,10])
   ax1.plot(U_metek_ts.values,U_gill_ts.values+U_residuals,'r+', label='U (m/s)')
   ax1.legend(loc=2)

   ax2 = ax1.twinx()
   pylab.ylim([-10,10])
   pylab.xlim([-10,10])
   ax2.plot(V_metek_ts.values,V_gill_ts.values+V_residuals,'b+', label='V (m/s)')
   ax2.plot(W_metek_ts.values,W_gill_ts.values+W_residuals,'g+', label='W (m/s)')
   ax2.plot([-10,10],[-10,10],'k-')
   ax2.legend(loc=1)
   pylab.savefig('sonic_intercomparison.png', dpi=300)



   sys.exit()
      
if __name__=='__main__':
   main()
   
   
   
