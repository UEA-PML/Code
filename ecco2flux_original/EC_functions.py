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
from numpy import copy, nanmean
import pylab 
import math
import argparse
import pandas as pd
from scipy import signal

# Function: pol2cart(phi, rho)
# Author: Taken from StackOverflow
# Description:
#  Converts from polar to cartesian co-ordinates
#
# Inputs:
#  phi - floating point degrees
#  rho - the "radius"  
# Outputs:
#  x and y co-ordinates
# Call function as:
#  x, y = pol2cart(phi, rho) 
def pol2cart(phi, rho):
    x = rho * math.cos(math.radians(phi))
    y = rho * math.sin(math.radians(phi))
    return(x, y)

# Function: cart2pol(x, y)
# Author: Taken from StackOverflow
# Description:
#  Converts from cartesian to polar co-ordinates
#
# Inputs:
#  x and y co-ordinates
# Outputs:
#  phi - floating point degrees
#  rho - the "radius"
# Call function as:
#  phi, rho = cart2pol(x, y) 
def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = math.degrees(np.arctan2(y, x))
    return(phi, rho)

# Function: angular_offset(U_raw, V_raw, W_raw, alpha, beta)
# Author: Tim Smyth
# Description: 
#  Function to correct wind data from Sonic 3D-wind anemometer for any angular offsets
#  between the motionpak and sonic
#
# Inputs:
#  u_ms - N element array of floating point horizontal vector wind (u) in m/s
#  v_ms - N element array of floating point horizontal vector wind (v) in m/s
#  w_ms - N element array of floating point vertical vector wind (w) in m/s
#  alpha - port / starboard orientation in radians 
#  beta - stern / bow orientation in radians
# Outputs:
#  u_ms - N element array of floating point horizontal vector wind (u) in m/s
#  v_ms - N element array of floating point horizontal vector wind (v) in m/s
#  w_ms - N element array of floating point vertical vector wind (w) in m/s
#
# Call function as:
# u_ms, v_ms, w_ms = angular_offset(u_ms, v_ms, w_ms, alpha, beta, polarity)  

def angular_offset(U_raw, V_raw, W_raw, alpha, beta, polarity):
   from numpy import cos, sin

   # reverse polarity of U axis for rotation below
   if (polarity == 1):
      U_raw = -U_raw
   
   u_ms = U_raw*cos(beta)*cos(alpha) + V_raw*sin(alpha)*cos(beta) + W_raw*sin(beta)
   v_ms = -U_raw*sin(alpha) + V_raw*cos(alpha)
   w_ms = -U_raw*sin(beta)*cos(alpha) - V_raw*sin(alpha)*sin(beta) + W_raw*cos(beta)

   #reverse U axis polarity again
   if (polarity == 1):
      u_ms = -u_ms

   return u_ms, v_ms, w_ms

"""
Functions for eddy correlation flux calculations.

Translated and vectorized from NOAA bulk model and flux processing MATLAB scripts.

Uncomment code at the end of this file and execute '%run flux.py' from the iPython command line to test functions.

Byron Blomquist, CU/CIRES, NOAA/ESRL/PSD3
v1: Oct 2014
"""

def motcorr(son, mot, heading, sens_disp, fsam = 10):
    """
    Motion corrections for ship-based flux measurements, after Edson 98.

    Usage: acc, uvw, uvwplat, xyzplat = motcorr(son,mot,heading,sens_disp)

    Inputs:
        son: Nx4 array of time and sonic wind velocities.
            N is typically 36000 for 10 Hz data.  See read_son().
        mot: Nx7 array of time, linear acceleration and angular rate variables.
            See read_mot().
        heading: N array of ship heading in radians--usually from a laser ring gyro.
        sens_disp: 3 element vector (x,y,z) of the measured displacement
            between the motion pak and sonic sensor volume: +x forward, +y port, +z up.
        fsam (optional): sampling rate, default is 10 Hz.

    Outputs:
		acc: Nx3 linear accelerations in Earth reference 
        uvw: Nx3 array of motion corrected wind components (u,v,w)
            where +u is wind velocity toward the bow and +v is the crosswind
            component from port (left side).  +w is the upward component.
        uvwplat: Nx3 array of platform velocity components for the
            same axes.
        xyzplat: Nx3 array of platform displacements on x,y,z axes, in meters.

    """

    import numpy as np
    from numpy import sqrt, sum, mean, unwrap, size, cross

    son_vel = son[:,1:4]    # sonic wind velocities
    plat_acc = mot[:,1:4]   # platform accelerations
    plat_rate = mot[:,4:]   # platform rotational rates

    # gravity, m/s2
    grav = sqrt(sum(mean(plat_acc,axis=0)**2))

    # unwrapped slow rotation about z axis
    heading = unwrap(heading)

    # euler angles
    euler = angles(plat_acc, plat_rate, heading, grav, sf=fsam)
    #np.savetxt('euler.txt',euler,fmt='%.18e',delimiter='\t')

    # platform velocities
    acc, uvwplat, xyzplat = accels(plat_acc, euler, grav, sf=fsam)

    # Calculate winds in earth based frame.  The angular velocity calculated as
    # the cross product between the angular rate vector and position vector.
    # The measured and angular velocities are in the platform frame and must be
    # rotated into the earth frame.
    #   uvw = measured velocity + angular rate induced velocities + integrated
    #       accelerometers.

    # sensor displacement array
    R = np.ones((size(uvwplat, 0), 3))
    R[:,0] *= sens_disp[0]
    R[:,1] *= sens_disp[1]
    R[:,2] *= sens_disp[2]

    # correct for angular velocities from sensor displacement
    uvw = son_vel + cross(plat_rate, R)

    # rotate to earth coordinates
    uvw = trans(uvw, euler)

    # Compute lag between W and platform vertical velocity.
    # This should be at most a few data points.
    lags = lagcorr(uvw[0:3000,2], uvwplat[0:3000,2], 20)

    # Now sum uvw with lag-shifted uvwplat to yield correct winds.
    # uvw are in earth coordinates, corrected for platform pitch,
    # roll, yaw and heave.  No correction for ship speed at this point.
    N = np.size(uvw, axis=0)    # number of points in timeseries
    if lags >= 0:    # plat lagged w/ respect to wind
        uvw[0:N-lags,:] = uvw[0:N-lags,:] + uvwplat[lags:N,:]
    elif lags < 0:  # wind lagged w/ respect to plat
        uvw[-lags:N,:] = uvw[-lags:N,:] + uvwplat[0:N+lags,:]

    return (acc, uvw, uvwplat, xyzplat)


def angles(accm, ratem, heading, grav, sf=10):
    """
    Computes euler angles from platform angular rates and accelerations.
    Uses the gravity vector to determine slow pitch and roll.

    The angles are estimated from:
        angle = slow_angle (from accelerometers) + fast_angle (integrated rate
            sensors)

    Usage: euler = angles(accm, ratem, heading, grav)

    Inputs:
        accm = Nx3 array of platform accelerations.  See read_mot().
        ratem = Nx3 array of platform rotation rates.  See read_mot().
        heading = N vector of unwrapped platform heading (slow psi angle) in
            radians.
        grav: gravity computed from accelerometers, m/s2.
        sf (optional) = sampling rate, default is 10 Hz.

    Output:
        euler = Nx3 array of euler angles

    """

    import numpy as np
    import scipy.signal as sig
    import scipy.integrate as int
    #*** import sig_proc as sp

    # Define high pass filter
    wp = 1/25./(sf/2.)
    ws = 0.8*wp
    n,wn = sig.buttord(wp, ws, 3, 7)
    bhi, ahi = sig.butter(n, wn, btype ='highpass')

    rates = sig.detrend(ratem, axis=0)  # detrend angular rates

    # slow roll
    phiraw = np.arctan2(accm[:,1], grav)
    phislow = phiraw - sig.filtfilt(bhi, ahi, phiraw, padtype='constant', padlen=1000)
    #*** phislow = phiraw - sp.filter_padded(phiraw, bhi, ahi) # pad with data

    # slow pitch
    thetaraw = np.arctan2(accm[:,0], grav)
    thetaslow = -thetaraw - sig.filtfilt(bhi, ahi, -thetaraw, padtype='constant', padlen=1000)
    #*** thetaslow = -thetaraw - sp.filter_padded(-thetaraw, bhi, ahi) # pad with data

    # slow yaw
    psislow = -heading - sig.filtfilt(bhi, ahi, -heading, padtype='constant', padlen=1000)   # slow yaw
    #*** psislow = -heading - sp.filter_padded(-heading, bhi, ahi) # pad with data

    #np.savetxt('phislow.txt',phislow,fmt='%.18e',delimiter='\t')
    #np.savetxt('thetaslow.txt',thetaslow,fmt='%.18e',delimiter='\t')
    #np.savetxt('psislow.txt',psislow,fmt='%.18e',delimiter='\t')

    # Integrate and filter angle rates, then add to slow angles
    #rate_int = np.cumsum(rates, axis=0) - 0.5*rates - 0.5*rates[0,:]/sf
    rate_int = int.cumtrapz(rates, dx=1./sf, axis=0, initial=0)
    phi = phislow + sig.filtfilt(bhi, ahi, rate_int[:,0], padtype='constant', padlen=1000)
    theta = thetaslow + sig.filtfilt(bhi, ahi, rate_int[:,1], padtype='constant', padlen=1000)
    psi = psislow + sig.filtfilt(bhi, ahi, rate_int[:,2], padtype='constant', padlen=1000)
    #*** phi = phislow + sp.filter_padded(rate_int[:,0], bhi, ahi)      # roll
    #*** theta = thetaslow + sp.filter_padded(rate_int[:,1], bhi, ahi)  # pitch
    #***psi = psislow + sp.filter_padded(rate_int[:,2], bhi, ahi)      # yaw

    euler = np.column_stack((phi,theta,psi))

    return euler


def accels(accm, euler, grav, sf=10):
    """
    Integrate linear accelerations to get platform velocity
        and displacement. After each integration, signals are
        high pass filtered to remove low frequency effects.

    Usage: acc, uvwplat, xyzplat = accels(accm, euler, grav)

    Inputs:
        accm = Nx3 array of platform accelerations.  See read_mot().
        euler = Nx3 array of euler angles.  See angles().
        grav: gravity computed from accelerometers, m/s2.
        sf (optional) = sampling frequency in Hz.  Default is 10Hz.

    Outputs:
		acc		= Nx3 linear accelerations in Earth reference 
        uvwplat = Nx3 linear velocities at the point of motion measurement.
        xyzplat = Nx3 platform displacements from mean position.

"""
    import numpy as np
    import scipy.signal as sig
    import scipy.integrate as int
    #*** import sig_proc as sp

    # Define high pass filter which retains real acceleration but removes drift
    wp = 1/25./(sf/2.)
    ws = 0.8*wp
    n,wn = sig.buttord(wp, ws, 3, 7)
    bhi, ahi = sig.butter(n, wn, btype ='highpass')

    acc = trans(accm, euler)    # first rotate into geophysical coordinates
    acc[:,2] = acc[:,2] - grav  # remove gravity from z axis accelerations
#     np.savetxt('acc.txt',acc,fmt='%.18e',delimiter='\t')

    # High-pass filter accelerations
    acc_filt = sig.filtfilt(bhi, ahi, acc, axis=0)
    #*** acc = sp.filter_padded(acc, bhi, ahi)
#     np.savetxt('acc_filt.txt',acc,fmt='%.18e',delimiter='\t')

    # Integrate accelerations to get velocity
#    uvwplat = np.cumsum(acc, axis=0) - 0.5*acc - 0.5*acc[0,:]/sf
    uvwplat = int.cumtrapz(acc_filt, dx=1./sf, axis=0, initial=0)

    # Filter velocities to remove drift
    uvwplat = sig.filtfilt(bhi, ahi, uvwplat, axis=0)
    #*** uvwplat = sp.filter_padded(uvwplat, bhi, ahi)

    # Integrate again to get displacements
#    xyzplat = np.cumsum(uvwplat, axis=0) - 0.5*uvwplat - 0.5*uvwplat[0,:]/sf
    xyzplat = int.cumtrapz(uvwplat, dx=1./sf, axis=0, initial=0)

    # Filter displacements to remove drift
    xyzplat = sig.filtfilt(bhi, ahi, xyzplat, axis=0)
    #*** xyzplat = sp.filter_padded(xyzplat, bhi, ahi)

    return (acc, uvwplat, xyzplat)


def dbl_rot(uvw):
    """
    Performs 2-angle rotation of uvw wind vectors into the mean flow
    streamlines such than mean v and mean w are zero.  The azimuth
    rotation here is clockwise to conform with the meteorological convention.

    Usage: uvw_str, tilt, azm = dbl_rot(uvw)

    Input: uvw = Nx3 numpy array of wind components.

    Output: uvw_str = Nx3 array of rotated wind vectors
            tilt = mean tilt angle in deg
            azm = mean relative wind direction or azimuth in deg
            ** Modification - TJS (11/12/2018)
            Umean = mean of U wind vector for averaging period
            Vmean = mean of V wind vector for averaging period
            Wmean = mean of W wind vector for averaging period

    """
    from numpy import sin, cos, arctan2, sqrt, column_stack, copy, pi, nanmean
    #from stats import nanmean

    u = copy(uvw[:,0])
    v = copy(uvw[:,1])
    w = copy(uvw[:,2])

    Ub = nanmean(u)
    Vb = nanmean(v)
    Wb = nanmean(w)
    Sb = sqrt(Ub**2 + Vb**2)

    tilt = arctan2(Wb, Sb)
    azm = arctan2(Vb, Ub)

    Ustr = +u*cos(azm)*cos(tilt) + v*sin(azm)*cos(tilt) + w*sin(tilt)
    Vstr = -u*sin(azm)      +      v*cos(azm)
    Wstr = -u*cos(azm)*sin(tilt) - v*sin(azm)*sin(tilt) + w*cos(tilt)

    Umean = nanmean(Ustr)
    Vmean = nanmean(Vstr)
    Wmean = nanmean(Wstr)

    uvw_str = column_stack((Ustr,Vstr,Wstr))

    return (uvw_str, tilt*180/pi, azm*180/pi, Umean, Vmean, Wmean)


def azm_rot(uvw, azm, ccw=True):
    """
    2-D coordinate rotation - rotates wind azimuth by the specified angle.

    Usage: uvw_rot = azm_rot(uvw, azm)

    Inputs: uvw = Nx3 numpy array of wind components.
            azm = rotation angle in radians
            ccw = boolean indicating direction of rotation.  Default is True for
                    ccw rotation, which is the normal mathematical convention.
                    False specifies a clockwise rotation (-alpha) which is the
                    meteorological convention for wind direction.  Setting ccw
                    True is the reverse azimuth rotation in double_rotation(uvw).

    Output: uvw_rot = Nx3 array of rotated wind vectors

    """
    from numpy import sin, cos, arctan2, sqrt, column_stack, copy, pi
    from stats import nanmean

    u = copy(uvw[:,0])
    v = copy(uvw[:,1])

    if ccw:
        Urot =  u*cos(azm) - v*sin(azm)
        Vrot =  u*sin(azm) + v*cos(azm)
    else:
        Urot =  u*cos(azm) + v*sin(azm)
        Vrot = -u*sin(azm) + v*cos(azm)

    uvw_rot = column_stack((Urot,Vrot,uvw[:,2]))

    return uvw_rot


def trans(src, angles, inv=False):
    """
    This subroutine rotates a vector from one cartesian basis to
    another, based upon the three Euler angles, defined as
    rotations around the reference axes, xyz. The axes in the
    rotated frame are x'y'z'.

    Usage: dest = trans(src, angles)

    Python version of June 27, 1996 Edson matlab code

    Input 'src' is a 3 col numpy array of vector components, x y z.
    Input 'angles' is a 3 col numpy array of phi, theta, psi, where:
        phi   = angles[:,0] - rotation of x'y'z' about x axis (roll)
        theta = angles[:,1] - rotation of x'y'z' about y axis (pitch)
        psi   = angles[:,2] - rotation of x'y'z' about z axis (yaw)

    Output 'dest' is the 3 column rotated vector components.

    The transformation is given by

    u' = [A]u      OR     u = [transpose(A)]u'

    Optional input 'inv' specifies the direction of transform:
    inv = False  "src" is transformed from x'y'z' to xyz (default)
    inv = True  "src" is transformed from xyz to x'y'z'

    For rotation of wind measurements from body coordinates to earth
    coordinates 'inv' is zero, u is along the x axis, v is y axis and
    w is z axis.

    This is a 321 rotation. This means that the first rotation
    is around the 3 axis (z-axis, angle psi), the second rotation is then
    about the intermediate 2 axis (y-axis, theta), and the third rotation
    is about the intermediate 1 axis (x-axis, phi).

"""
    from numpy import size, sin, cos, column_stack

    if size(src,0) == 3:
        src = src.T # transpose from 3 row to 3 column array
    if size(angles,0) == 3:
        angles = angles.T # ditto

    p  = angles[:,0]
    t  = angles[:,1]
    ps = angles[:,2]

    up = src[:,0]
    vp = src[:,1]
    wp = src[:,2]

    if inv:
        u = up*cos(t)*cos(ps) + vp*cos(t)*sin(ps) - wp*sin(t)
        v = (up*(sin(p)*sin(t)*cos(ps) - cos(p)*sin(ps)) +
            vp*(sin(p)*sin(t)*sin(ps) + cos(p)*cos(ps)) + wp*(cos(t)*sin(p)))
        w = (up*(cos(p)*sin(t)*cos(ps) + sin(p)*sin(ps)) +
            vp*(cos(p)*sin(t)*sin(ps) - sin(p)*cos(ps)) + wp*(cos(t)*cos(p)))
    else:
        u = (up*cos(t)*cos(ps) + vp*(sin(p)*sin(t)*cos(ps) -
            cos(p)*sin(ps)) + wp*(cos(p)*sin(t)*cos(ps) + sin(p)*sin(ps)))
        v = (up*cos(t)*sin(ps) + vp*(sin(p)*sin(t)*sin(ps) +
            cos(p)*cos(ps)) + wp*(cos(p)*sin(t)*sin(ps) - sin(p)*cos(ps)))
        w = up*(-sin(t)) + vp*(cos(t)*sin(p)) + wp*(cos(t)*cos(p))

    dest = column_stack((u,v,w))
    return dest


def lagcorr(X, Y, maxlags):
    """
    Computes time lag between two signals via cross-correlation.
    Assumes correlation between signals is sufficiently strong.
    Uncorrelated signals will not yield a meaningful result.

    Usage: lags = lagcorr(X, Y, maxlags)

    Inputs:
        X, Y = input 1-D arrays of equal length.
        maxlags = correlation is computed over range of +/- maxlags data points.

    Output:
        lags = number of data points Y is lagged with respect to X.
            Negative value means X is lagged with respect to Y.

    """
    import numpy as np
    import scipy.signal as sig

    N = np.size(X)
    xycorr = sig.correlate(X, Y, mode='same')
    # select segment centered on zero lag
    midIdx = np.size(xycorr)/2  # midpoint index == zero lag
    seg = np.copy(xycorr[midIdx-maxlags : midIdx+maxlags])

    minIdx = seg.argmin()  # index of minimum correlation
    maxIdx = seg.argmax()  # index of maximum correlation
    # use correlation with greatest magnitude (pos or neg)
    if np.abs(seg.max()) > abs(seg.min()):
        lags = maxIdx - maxlags
    else:
        lags = minIdx - maxlags

    return lags


def compute_q_T(Ts, W_lic, T_bar, Q_bar, W_bar):
    """
    computes high rate q' and T' from raw sonic Ts and Licor H2O mmol/m3

    usage: q_fast, T_fast = compute_q_T(Ts, W_lic, T_bar, Q_bar, W_bar)

    inputs: Ts: 10 Hz sonic temperature, C
            W_lic: 10 Hz Licor H2O number density, mmol/m3
            T_bar: mean T from T/RH sensor, C
            Q_bar: mean specific humidity from T/RH sensor, g/kg
            W_bar: mean H2O number density, mmol/m3

    note: the number of elements in W-Bar and T_bar determines the number of
            fast data segments. i.e. W_bar.shape[0]=6 divides fast data into
            6 segments (i.e. 10-min segments for hourly input data)
            W_bar and T_bar must be the same size.

    outputs: q_fast: 10 Hz specific humidity, g/kg
             T_fast: 10 Hz temperature, C

    """
    import numpy as np

    if T_bar.shape[0] != W_bar.shape[0]:
        raise ValueError, 'compute_q_T: W_bar.shape[0] != T_bar.shape[0]'
    if Ts.shape[0] != W_lic.shape[0]:
        raise ValueError, 'compute_q_T: Ts.shape[0] != W_lic.shape[0]'

    N = T_bar.shape[0]                  # number of segments
    Npts = Ts.shape[0]/N                # points per segment
    if Npts*N != Ts.shape[0]:
        print 'Warning: compute_q_T: fast array not evenly divisible'

    T_bar_K = np.copy(T_bar) + 273.15   # convert to deg K
    Ts_K = np.copy(Ts) + 273.15         # convert to deg K
    q_bar = Q_bar/1000.0                # convert to kg/kg
    w_bar = W_bar/1000.0                # mean H2O number density, moles/m3
    w_lic = np.copy(W_lic)/1000.0       # convert to moles/m3
    T_fast = np.zeros_like(Ts)*np.nan   # empty output arrays
    q_fast = np.zeros_like(W_lic)*np.nan

    for ii,jj in zip(np.arange(0,Ts.shape[0],Npts),np.arange(N)):
        # solving two equations with two unknowns via inverse of 2x2 matrix
        # (w'/w) = a*q' + b*T'
        # Ts' = c*q' + d*T'
        # where A=[[a,b],[c,d]], defined below
        a = 1/(q_bar[jj]*(1 + (0.378/0.622)*q_bar[jj]))
        b = -1/T_bar_K[jj]
        c = 0.51*T_bar_K[jj]
        d = 1 + 0.51*q_bar[jj]
        det = a*d - b*c     # matrix determinant
        # solve [x1,x2] = [A]^-1 [y1,y2]
        q_fast[ii:ii+Npts] = (d/det)*(w_lic[ii:ii+Npts]/w_bar[jj]) + (-b/det)*Ts_K[ii:ii+Npts]
        T_fast[ii:ii+Npts] = (-c/det)*(w_lic[ii:ii+Npts]/w_bar[jj]) + (a/det)*Ts_K[ii:ii+Npts]

    q_fast *= 1000.0    # back to g/kg units
    T_fast -= 273.15
    return (q_fast,T_fast)


def psi_fu_c(zet):
    """
    usage psi = psi_fu_c(zet)

    input: zet: z/L Monin-Obukhov stability parameter (unitless)

    """
    import numpy as np
    import util

    psi = 3.9*(1.0/(1 - 20*zet)**0.333 - zet - 1.0/(7.0 - zet))**(2.0/3.0)
    uns = util.find(zet>0)
    if uns.size > 0:
        psi[uns] = 0.92*3.9*(1.0 + 2.5*zet[uns]**(2.0/3.0))

    return psi


def psi_ft_c(zet):
    """
    usage psi = psi_ft_c(zet)

    input: zet: z/L Monin-Obukhov stability parameter (unitless)

    """
    import numpy as np
    import util

    psi = 5.5/(1.0 - 7.0*zet)**(2.0/3.0)
    uns = util.find(zet>0)
    if uns.size > 0:
        psi[uns] = 5.5*(1.0 + 2.5*zet[uns]**(2.0/3.0))

    return psi


def idiss_struc_func(Pxx, Fxx, fmin, fmax, rwspd):
    """
    Compute inertial dissipation structure function parameter from median

    usage Cx2 = idiss_struc_func(Pxx, Fxx, fmin, fmax, rwspd)

    inputs: Pxx: variance spectrum, smoothed / bin averaged
                spectrum may be single vector or 2-D with one spectrum per row
            Fxx: frequency axis for spectrum, 1-D same length as Pxx
            fmin: start frequency
            fmax: end frequency
            rwspd: relative wind speed, vector with one value for each
            spectrum in Pxx

    output: Cx2: structure function parameter as nparray

    """
    import numpy as np
    import util

    if Pxx.ndim > 1:
        if Pxx.shape[0] != rwspd.shape[0]:
            raise ValueError, 'idiss_struc_func: rwspd.shape[0] != Pxx.shape[0]'

    cS = 4.0
    Cu2 = np.zeros(Pxx.shape[0])
    fmax += 1 # so range indexing includes fmax
    for ii in np.arange(6):
        Cu2[ii] = np.median(cS*((2.*np.pi/rwspd[ii])**(2.0/3.0)*Pxx[ii,fmin:fmax]*
                Fxx[fmin:fmax]**(5.0/3.0)))

    return Cu2


def idiss_spec_fit(Pxx, Fxx, fmin, fmax, rwspd):
    """
    Compute the intertial dissipation structure function parameter with a
        spectrum fit.

    usage Cx2 = idiss_spec_fit(Pxx, Fxx, fmin, fmax, rwspd)

    inputs: Pxx: variance spectrum, smoothed / bin averaged
                spectrum may be single vector or 2-D with one spectrum per row
            Fxx: frequency axis for spectrum, 1-D same length as Pxx
            fmin: start frequency for fit
            fmax: end frequency for fit
            rwspd: relative wind speed, vector with one value for each spectrum in Pxx

    output: Cx2: structure function parameter as ndarray
    """
    import numpy as np
    import util

    if Pxx.ndim > 1:
        if Pxx.shape[0] != rwspd.shape[0]:
            raise ValueError, 'idiss_spec_fit: rwspd.shape[0] != Pxx.shape[0]'

    cS = 4.0
    Cx2 = np.zeros(Pxx.shape[0])
    fmax += 1 # so range indexing includes fmax
    # fit spectrum
    logf = np.log(Fxx[fmin:fmax])
    logP = np.zeros((Pxx.shape[0],logf.shape[0]))
    for ii in np.arange(logP.shape[0]): # compute log of spectra
        logP[ii,:] = np.log(cS*((2*np.pi/rwspd[ii])**(2.0/3.0)*
                Pxx[ii,fmin:fmax]*Fxx[fmin:fmax]**(5.0/3.0)))
    zpoly = np.polyfit(logf, logP.T, 4)
    zslope = np.copy(zpoly[3,:])

    ii = util.find(zslope<0.25)
    Cx2[ii] = np.exp(zpoly[4,ii])

    jj = util.find(zslope>=0.25)
    zz = np.copy(jj)
    plus = np.zeros(zslope.shape[0], dtype=int)
    if jj.shape[0] > 0:
        while 1: # iteratively increase fmin and recompute slopes
        #Zslope = 4*Zpoly(1).*(log(Fson_uj(fminU+plus)).^3)+3*Zpoly(2).*
        #   (log(Fson_uj(fminU+plus)).^2)+2*Zpoly(3).*log(Fson_uj(fminU+plus))+
        #   Zpoly(4)
            zold = np.copy(zslope)
            zslope[jj] = (4.0*zpoly[0,jj]*(np.log(Fxx[fmin + plus[jj]])**3) +
                          3.0*zpoly[1,jj]*(np.log(Fxx[fmin + plus[jj]])**2) +
                          2.0*zpoly[2,jj]*(np.log(Fxx[fmin + plus[jj]])) +
                              zpoly[3,jj])
            plus[jj] += 1
            jj = util.find((zslope <= zold)*(np.abs(zslope) >= 0.2)*
                    (plus+fmin < fmax)) # find slopes where all conditions are True
            if jj.shape[0] == 0:        # break loop when none meet the above criteria
                break
        Cx2[zz] = np.exp(zpoly[0,zz]*(np.log(Fxx[fmin+plus[zz]])**4) +
                         zpoly[1,zz]*(np.log(Fxx[fmin+plus[zz]])**3) +
                         zpoly[2,zz]*(np.log(Fxx[fmin+plus[zz]])**2) +
                         zpoly[3,zz]*np.log(Fxx[fmin+plus[zz]]) +
                         zpoly[4,zz])
    return Cx2


def idiss(Cu2, Ct2, Cq2, qbar, Tbar, zu, zq, L_bulk, tsr_bulk, qsr_bulk):
    """
    Inertial dissipation calculations

    usage: usid, tsid, qsid, Lid = idiss(Cu2, Ct2, Cq2, qbar, Tbar, zu, zq,
                                            L_bulk, tsr_bulk, qsr_bulk)

    inputs may be single values or vectors of equal length
    inputs: Cu2, Ct2, Cq2: structure function parameters
            qbar, Tbar: mean q and T from bulk T/RH sensors
            zu, zq: wind and q measurement heights (m)
            L_bulk, tsr_bulk, qsr_bulk: bulk model results

    outputs: usid, tsid, qsid, Lid

    """
    import numpy as np

    # set up constants
    vk = 0.4                        # von Karman const
    alpahu = 0.52
    c_coeff = 4*vk/alpahu**0.667
    g = 9.82                        # gravity
    Cp = 1004                       # heat capacity dry air
    Cle = 2500 - 2.274*Tbar         # latent heat evap
    T = Tbar + 273.15               # deg K
    Tv = T*(1 + 0.61*(qbar/(1-qbar))) # virtual T
    rho = 1013.25/287/Tv            # air density
    zet_u = zu/L_bulk               # z/L for u and q
    zet_q = zq/L_bulk

    usid = np.sqrt(Cu2*zu**(2.0/3.0)/psi_fu_c(zet_u))
    tsid = np.sign(zet_u)*np.sqrt(Ct2*zu**(2.0/3.0)/psi_ft_c(zet_u))
    qsid = np.sqrt(Cq2*zq**(2.0/3.0)/psi_ft_c(zet_q))
    tsid = tsid - 0.00051*T*qsid/rho

    us0 = usid
    ts0 = tsr_bulk
    qs0 = qsr_bulk

    for ii in np.arange(10):
        zet = vk*zu*g/Tv*(ts0 + 0.00061*T*qs0/rho)/us0**2
        us0 = np.sqrt(Cu2*zu**(2.0/3.0)/psi_fu_c(zet))

    Lid = zu/zet

    return (us0, tsid, qsid, Lid)

def reject_outliers(data,m):
    stdev = np.std(data)
    mean = np.mean(data)
    maskMin = mean - stdev * m
    maskMax = mean + stdev * m
    test = np.where(data<=maskMin) or np.where(data>=maskMax)
    test = test[0]
    if (len(test) > 0):
       for i in range(len(test)):
          data[test[i]] = np.nan
    return data
    
def reject_ones(data):
    last_switch_index = 0
    residual_gas = 0
    
    test = np.where(data != 0)
    test = test[0]
    if (len(test) > 0):
       for i in range(len(test)):
          switch_flag = 1
          last_switch_index = test[i]
          data[test[i]] = np.nan
       residual_gas = last_switch_index + 100
       if (residual_gas > len(data)):
          residual_gas = len(data)
       for j in range(last_switch_index, residual_gas):
          data[j] = np.nan
          
    return data
    
def platmot_corrected_gas_flux(input_gas, valve_switch_state, Zplat, Wplat, Aplat, Wprime):

   # linear fit of gas vs time
   # Create an array of same size as input_gas which increments at 0.1
   time_s_array = np.zeros(len(input_gas))*0.
   for ss in range(len(input_gas)):
      time_s_array[ss] = float(ss)*0.1
   
   # reject outliers from the gas concentration datasets
   input_gas = reject_outliers(input_gas.values,5)

   # reject outliers from the platform movement datasets
   Zplat = reject_outliers(Zplat.values,5)
   Wplat = reject_outliers(Wplat.values,5)
   Aplat = reject_outliers(Aplat.values,5)
   
   # reject values when the valve is triggered
   valve_switch_state = reject_ones(valve_switch_state)

   # identify the elements in the matrix which are NaN and therefore an outlier
   idy = np.where(np.isfinite(input_gas) & np.isfinite(Zplat) & np.isfinite(Wplat) & np.isfinite(Aplat) & np.isfinite(valve_switch_state))
   
   # calculate statistics to return to main routine on gas concentration
   mean_input_gas_conc = np.mean(input_gas[idy])
   std_input_gas_conc = np.std(input_gas[idy])
   
   # rather than doing a simple mean, produce a linear fit
   m_Zplat,Zplat_c = np.polyfit(time_s_array[idy],Zplat[idy],1)
   m_Wplat,Wplat_c = np.polyfit(time_s_array[idy],Wplat[idy],1)
   m_Aplat,Aplat_c = np.polyfit(time_s_array[idy],Aplat[idy],1)

   Zplat_lin = Zplat_c + m_Zplat*time_s_array[idy]
   Wplat_lin = Wplat_c + m_Wplat*time_s_array[idy]
   Aplat_lin = Aplat_c + m_Aplat*time_s_array[idy]

   Wprime_gas = Wprime[idy]

   Zplat_prime = Zplat[idy] - Zplat_lin
   Wplat_prime = Wplat[idy] - Wplat_lin
   Aplat_prime = Aplat[idy] - Aplat_lin
   
   # against displacement
   m_gas,gas_c = np.polyfit(time_s_array[idy],input_gas[idy],1)
   input_gas_lin = gas_c + m_gas*time_s_array[idy]
   input_gas_prime = input_gas[idy] - input_gas_lin

   Zplat_cov = input_gas_prime*Zplat_prime
   mu_Zplat = np.mean(Zplat_cov)/np.var(Zplat_prime)
   input_gas = input_gas[idy] - mu_Zplat*Zplat_prime

   # regenerate the time matrix based on the number of accepted elements
   time_s_array = np.zeros(len(input_gas))*0.
   for ss in range(len(input_gas)):
      time_s_array[ss] = float(ss)*0.1

   # against velocity
   # recompute c'
   m_gas,gas_c = np.polyfit(time_s_array,input_gas,1)
   input_gas_lin = gas_c + m_gas*time_s_array
   input_gas_prime = input_gas - input_gas_lin

   Wplat_cov = input_gas_prime*Wplat_prime
   mu_Wplat = np.mean(Wplat_cov)/np.var(Wplat_prime)
   input_gas = input_gas - mu_Wplat*Wplat_prime

   # against acceleration
   # recompute c'   
   m_gas,gas_c = np.polyfit(time_s_array,input_gas,1)
   input_gas_lin = gas_c + m_gas*time_s_array
   input_gas_prime = input_gas - input_gas_lin
   
   Aplat_cov = input_gas_prime*Aplat_prime
   mu_Aplat = np.mean(Aplat_cov)/np.var(Aplat_prime)
   input_gas = input_gas - mu_Aplat*Aplat_prime
   
   m_gas,gas_c = np.polyfit(time_s_array,input_gas,1)
   input_gas_lin = gas_c + m_gas*time_s_array
   input_gas_prime = input_gas - input_gas_lin

   input_gasF = np.mean(input_gas_prime*Wprime_gas)
   
   return mean_input_gas_conc, std_input_gas_conc, input_gasF

# Estimation of the longwave flux   
def longwave_flux(TairK,SST,Fc,P_licor_kPa,H2O_ppm):
   emiss_lw = 0.985       # Emissivity of the Ocean (Dickey et al.)
   sigmaSB = 5.6697e-8    # Stefan-Boltzmann 
   
   SSTK = SST + 273.15
   P_licor_hPa = 10.*P_licor_kPa
   
   ea = (P_licor_hPa*(H2O_ppm/1.0e+6))/(1.0 + (H2O_ppm/1.0e+6))
   
   # Berliand calculation of longwave flux
   lwest = -emiss_lw*sigmaSB*(TairK**4)*(0.39-0.05*np.sqrt(ea))*Fc - 4.0*emiss_lw*sigmaSB*(TairK**3)*(SSTK-TairK)

   # The convention for air-sea exchange is positive upwards
   return -lwest

# Determination of the mixing ratio (in g/kg)
def mixing_ratio(H2O_ppm,P_licor_kPa):
   B = 0.6219907 #kg/kg
   Pw = (P_licor_kPa*(H2O_ppm/1e+6))/(1.0 + (H2O_ppm/1e+6)) 
   X = B*Pw/(P_licor_kPa - Pw)
   
   return X
   
def alternative_Wprime(W):
   time_s_array = np.zeros(len(W))*0.
   for ss in range(len(W)):
      time_s_array[ss] = float(ss)*0.1
   
   m_W,W_c = np.polyfit(time_s_array,W,1)
   W_lin = W_c + m_W*time_s_array
   Wprime = W - W_lin

   return Wprime
   




