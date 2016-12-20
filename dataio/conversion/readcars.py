"""
Tools for parsing the CSIRO Atlas of Regional Seas (CARS) data

M. Rayson
University of Western Australia
Dec 2016
"""

import numpy as np 
import xray
from scipy.interpolate import interp2d, interp1d

from soda.utils.othertime import datetime64todatetime, YearDay

def load_cars_temp(carsfile, X, Y, Z, T):
    """
    Load data from the CARS
    """
    try:
        nx = X.shape[0]
        ny = Y.shape[0]
    except:
        raise Exception, 'X/Y need to be numpy array objects'

    assert nx == ny, 'X and Y should be vectors of the same length'

    cars = xray.open_dataset(carsfile)

    # Load the CARS coordinates
    Xcars = cars['mean'].lon.values
    Ycars = cars['mean'].lat.values

    # Convert to negative down
    Zmean = -1*cars['depth'].values
    Z_A = -1*cars['depth_ann'].values
    Z_SA = -1*cars['depth_semiann'].values

    nzmean = Zmean.shape[0]
    nzA = Z_A.shape[0]
    nzSA = Z_SA.shape[0]

    ## Calculate the annual and semi-annual harmonics
    #T_mean = cars['mean'][...]
    #A_SA = cars['an_cos'][...] + 1j*cars['an_sin'][...]
    #A_SSA = cars['sa_cos'][...] + 1j*cars['sa_sin'][...]

    ######
    # Interpolate the amplitude functions spatially

    # Construct output arrays
    Tmean_i = np.zeros((nzmean, ny))
    T_SAc_i = np.zeros((nzSA, ny))
    T_SAs_i = np.zeros((nzSA, ny))
    T_Ac_i = np.zeros((nzA, ny))
    T_As_i = np.zeros((nzA, ny))

    # Mean
    for kk in range(nzmean):
        cdata = cars['mean'][kk,...].values
        cdata[np.isnan(cdata)] = 0.
        Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        Tmean_i[kk,...] = Finterp(X,Y)

    # Annual
    for kk in range(nzA):
        cdata = cars['an_cos'][kk,...].values
        cdata[np.isnan(cdata)] = 0.
        Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        T_Ac_i[kk,...] = Finterp(X,Y)

        cdata = cars['an_sin'][kk,...].values
        cdata[np.isnan(cdata)] = 0.
        Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        T_As_i[kk,...] = Finterp(X,Y)

    # Annual
    for kk in range(nzSA):
        cdata = cars['sa_cos'][kk,...].values
        cdata[np.isnan(cdata)] = 0.
        Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        T_SAc_i[kk,...] = Finterp(X,Y)

        cdata = cars['sa_sin'][kk,...].values
        cdata[np.isnan(cdata)] = 0.
        Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        T_SAs_i[kk,...] = Finterp(X,Y)

    ######
    # Interpolate vertically
    kind='slinear'
    F = interp1d(Zmean, Tmean_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    Tmean_iz = F(Z)

    F = interp1d(Z_A, T_Ac_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    T_Ac_iz = F(Z)
    F = interp1d(Z_A, T_As_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    T_As_iz = F(Z)

    F = interp1d(Z_SA, T_SAc_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    T_SAc_iz = F(Z)
    F = interp1d(Z_SA, T_SAs_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    T_SAs_iz = F(Z)


    #######
    # Compute the temporal values from the amplitude
    tdt = datetime64todatetime(T)
    t_yday = YearDay(tdt)
    t_osc = 2*np.pi* t_yday / 366.


    # Example from: http://www.marine.csiro.au/~dunn/cars2009/
    #  Evaluate at day-of-year 45 (mid February)
    #  t = 2pi x 45/366 
    #  feb = mean + an_cos*cos(t) + an_sin*sin(t) + sa_cos*cos(2*t) + sa_sin*sin(2*t) 

    Tcars = Tmean_iz + T_Ac_iz*np.cos(t_osc[np.newaxis,...]) + \
        T_As_iz*np.sin(t_osc[np.newaxis,...]) + \
        T_SAc_iz*np.cos(2*t_osc[np.newaxis,...]) + \
        T_SAs_iz*np.sin(2*t_osc[np.newaxis,...])

    return Tcars


