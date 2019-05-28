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
from soda.utils.interpXYZ import interpXYZ

import pdb

def load_cars_temp(carsfile, X, Y, Z, T):
    """
    Load data from the CARS
    """
    try:
        nx = X.shape[0]
        ny = Y.shape[0]
    except:
        raise Exception('X/Y need to be numpy array objects')

    assert nx == ny, 'X and Y should be vectors of the same length'
    xlims = slice(np.min(X)-2, np.max(X)+2)
    ylims = slice(np.min(Y)-2, np.max(Y)+2)
    #indices = dict(lon=xlims,lat=ylims)

    cars = xray.open_dataset(carsfile)

    # Load the CARS coordinates
    xcars = cars['mean'].sel(lon=xlims,lat=ylims).lon.values.astype(np.float64)
    ycars = cars['mean'].sel(lon=xlims,lat=ylims).lat.values.astype(np.float64)
    Xcars,Ycars = np.meshgrid(xcars,ycars)

    # Convert to negative down
    Zmean = -1*cars['depth'].values.astype(np.float64)
    Z_A = -1*cars['depth_ann'].values.astype(np.float64)
    Z_SA = -1*cars['depth_semiann'].values.astype(np.float64)

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

    XYout = np.array([X,Y]).T
    #args = {'method':'kriging', 'NNear':8,'vrange':2.0}
    args = {'method':'idw', 'NNear':4,'p':2.0}
    #args = {'method':'curvmin'}
    # Mean
    for kk in range(nzmean):
        #print '%d of %d...'%(kk,nzmean)
        cmean = cars['mean'][kk,...]
        cdata = cmean.sel(lon=xlims,lat=ylims).values.squeeze()
        #cdata[np.isnan(cdata)] = 0.
        #Finterp = interp2d(Xcars, Ycars, cdata, kind='cubic')
        #Tmean_i[kk,...] = Finterp(X,Y)
        idx = ~np.isnan(cdata)
        if np.any(idx):
            Finterp = interpXYZ(np.array([Xcars[idx],Ycars[idx]]).T, XYout, **args)
            Tmean_i[kk,...] = Finterp(cdata[idx])

    # Annual
    for kk in range(nzA):
        #print '%d of %d...'%(kk,nzA)
        cdatax = cars['an_cos'][kk,...]
        cdata = cdatax.sel(lon=xlims,lat=ylims).values.squeeze()

        #cdata[np.isnan(cdata)] = 0.
        #Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        #T_Ac_i[kk,...] = Finterp(X,Y)

        idx = ~np.isnan(cdata)
        if np.any(idx):
            Finterp = interpXYZ(np.array([Xcars[idx],Ycars[idx]]).T, XYout, **args)
            T_Ac_i[kk,...] = Finterp(cdata[idx])

            cdatax = cars['an_sin'][kk,...]
            cdata = cdatax.sel(lon=xlims,lat=ylims).values.squeeze()
            T_As_i[kk,...] = Finterp(cdata[idx])

        #cdata[np.isnan(cdata)] = 0.
        #Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        #T_As_i[kk,...] = Finterp(X,Y)

    # Semi-Annual
    for kk in range(nzSA):
        #print '%d of %d...'%(kk,nzSA)
        ##cdata = cars['sa_cos'][kk,...].values
        ##cdata[np.isnan(cdata)] = 0.
        ##Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        ##T_SAc_i[kk,...] = Finterp(X,Y)

        ##cdata = cars['sa_sin'][kk,...].values
        ##cdata[np.isnan(cdata)] = 0.
        ##Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        ##T_SAs_i[kk,...] = Finterp(X,Y)

        #cdata = cars['sa_cos'][kk,...].values.squeeze()
        ##cdata[np.isnan(cdata)] = 0.
        ##Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        ##T_Ac_i[kk,...] = Finterp(X,Y)

        #idx = ~np.isnan(cdata)
        #Finterp = interpXYZ(np.array([Xcars[idx],Ycars[idx]]).T, XYout)
        #T_SAc_i[kk,...] = Finterp(cdata[idx])

        #cdata = cars['sa_sin'][kk,...].values.squeeze()
        #T_SAs_i[kk,...] = Finterp(cdata[idx])
        ##cdata[np.isnan(cdata)] = 0.
        ##Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        ##T_As_i[kk,...] = Finterp(X,Y)

        cdatax = cars['sa_cos'][kk,...]
        cdata = cdatax.sel(lon=xlims,lat=ylims).values.squeeze()

        #cdata[np.isnan(cdata)] = 0.
        #Finterp = interp2d(Xcars, Ycars, cdata, kind='linear')
        #T_Ac_i[kk,...] = Finterp(X,Y)

        idx = ~np.isnan(cdata)
        if np.any(idx):
            Finterp = interpXYZ(np.array([Xcars[idx],Ycars[idx]]).T, XYout, **args)
            T_SAc_i[kk,...] = Finterp(cdata[idx])

            cdatax = cars['sa_sin'][kk,...]
            cdata = cdatax.sel(lon=xlims,lat=ylims).values.squeeze()
            T_SAs_i[kk,...] = Finterp(cdata[idx])


    ######
    # Interpolate vertically
    def interp_z(Z, T, Zout, kind=2):
        F = interp1d(Z, T, kind=kind, axis=0,\
                bounds_error=False, fill_value='extrapolate')
        return F(Zout)

    Tmean_iz = interp_z(Zmean, Tmean_i, Z)
    T_Ac_iz = interp_z(Z_A, T_Ac_i, Z)
    T_As_iz = interp_z(Z_A, T_As_i, Z)
    T_SAc_iz = interp_z(Z_SA, T_SAc_i, Z)
    T_SAs_iz = interp_z(Z_SA, T_SAs_i, Z)
    
    #kind='slinear'
    #F = interp1d(Zmean, Tmean_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    #Tmean_iz = F(Z)

    #F = interp1d(Z_A, T_Ac_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    #T_Ac_iz = F(Z)
    #F = interp1d(Z_A, T_As_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    #T_As_iz = F(Z)

    #F = interp1d(Z_SA, T_SAc_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    #T_SAc_iz = F(Z)
    #F = interp1d(Z_SA, T_SAs_i, kind=kind, axis=0, bounds_error=False, fill_value=0.)
    #T_SAs_iz = F(Z)


    #######
    # Compute the temporal values from the amplitude
    try:
        tdt = datetime64todatetime(T)
    except:
        tdt = T # already datetime object

    t_yday = YearDay(tdt)
    t_osc = 2*np.pi* t_yday / 366.

    # Example from: http://www.marine.csiro.au/~dunn/cars2009/
    #  Evaluate at day-of-year 45 (mid February)
    #  t = 2pi x 45/366 
    #  feb = mean + an_cos*cos(t) + an_sin*sin(t) + sa_cos*cos(2*t) + sa_sin*sin(2*t) 

    Tcars = Tmean_iz[...,np.newaxis] + \
        T_Ac_iz[...,np.newaxis]*np.cos(t_osc[np.newaxis,...]) + \
        T_As_iz[...,np.newaxis]*np.sin(t_osc[np.newaxis,...]) + \
        T_SAc_iz[...,np.newaxis]*np.cos(2*t_osc[np.newaxis,...]) + \
        T_SAs_iz[...,np.newaxis]*np.sin(2*t_osc[np.newaxis,...])

    return Tcars


