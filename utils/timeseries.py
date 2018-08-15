# -*- coding: utf-8 -*-
"""
Collection of tools for plotting and analysis of time series data

Created on Wed Feb 06 12:37:17 2013
Author: Matt Rayson
Stanford University
"""

import numpy as np
import numexpr as ne
import matplotlib.pyplot as plt
from matplotlib import mlab
from matplotlib.lines import Line2D
from scipy import signal, interpolate
from datetime import datetime, timedelta
import operator
import xray
import pandas as pd

from . import othertime
from .uspectra import uspectra, getTideFreq
from .otherplot import stackplot
from .harmonic_analysis import harmonic_fit, harmonic_signal
from .mysignal import filt_gaps
from soda.dataio.netcdfio import queryNC

import pdb

class timeseries(object):
    """
    Class for handling time series data
    
    Methods include:
        - Power spectral density (plot)
        - Filtering
        - Interpolation
        - Plotting
    """
    
    basetime = datetime(1900,1,1)
    VERBOSE=False
    
    units=''
    long_name=''
    StationID = ''
    StationName = ''
    varname = ''
    Z=0.0
    X=0.0
    Y=0.0
    
    def __init__(self,t,y,**kwargs):
        
        self.__dict__.update(kwargs)        
        self.t = self._set_time(t) # independent variable (t,x, etc)
        self.y = y # dependent variable
        
        self.shape = self.y.shape
        self.ndim = len(self.shape)
        

        # Convert the time to seconds
        self.tsec = self._get_tsec(self.t)
        #if isinstance(self.t[0], np.datetime64):
        #    time = self.t
        #    self.tsec = ((time - time[0])*1e-9).astype(np.float64)
        #else:
        #    self.tsec = othertime.SecondsSince(self.t,basetime=self.basetime)
        
        self.ny = np.size(self.y)

        # make sure the original data is a masked array
        if not isinstance(self.y,np.ma.core.MaskedArray):
            mask = ~np.isfinite(self.y)
            self.y = np.ma.MaskedArray(self.y,mask=mask)
        
        self.dt, self.isequal = self._check_dt(self.tsec)

        self.Nt = self.t.shape[0]
        
        # Make sure that time is the last dimension
        if self.y.shape[-1] != self.Nt:
            self.y=self.y.T
        
    def psd(self, plot=True, nbandavg=1, scale=1.,**kwargs):
        """
        Power spectral density
        
        nbandavg = Number of fft bnns to average
        scale = scale factor (for plotting only)
        """
        
        if self.isequal==False and self.VERBOSE:
            print('Warning - time series is unequally spaced consider using lomb-scargle fft')
        

        NFFT = int(2**(np.floor(np.log2(self.ny/nbandavg)))) # Nearest power of 2 to length of data
        # Remove the mean and nan's
        y = self.y - self.y.mean()
        y[y.mask] = 0.

        Pyy,frq = mlab.psd(y,\
                Fs=2*np.pi/self.dt,\
                NFFT=NFFT,\
                window=mlab.window_hanning,\
                scale_by_freq=True)
        
        if plot:
            plt.loglog(frq,Pyy*scale,**kwargs)
            plt.xlabel('Freq. [$radians s^{-1}$]')
            plt.ylabel('PSD')
            plt.grid(b=1)
        
        return Pyy, frq

    def autocorr(self,normalize=False,axis=-1):
        """
        Autocorrelation calculation
        """

        assert self.isequal,\
             'Data must be equally spaced to perform this function.'

        N = self.ny
        M = int(N)/10

        ymean = self.y.mean(axis=axis)
        y = self.y - ymean
        k = list(range(1,M))
        tau = np.asarray(k,dtype=np.float)*self.dt

        Cyy = [1./(N-kk) * np.sum(y[...,0:-kk]*y[...,kk::],axis=axis) for kk in k ]

        if normalize:
            return Cyy/y.var(), tau
        else:
            return Cyy ,tau
            
        
    def specgram(self, NFFT=256,noverlap=128,plot=True,vv=29, **kwargs):
        """
        Spectrogram plot
        """
        from matplotlib.colors import LogNorm
        
        Pyy,frq,tmid = mlab.specgram(self.y-self.y.mean(),Fs=2*np.pi/self.dt,window=mlab.window_hanning,\
            scale_by_freq=True,noverlap=noverlap,NFFT=NFFT)
            
        if plot==True:
            ax3=plt.gca()
            plt.contourf(tmid*2*np.pi,frq,Pyy,vv,norm=LogNorm(),**kwargs)
            ax3.set_yscale('log')
            #ax3.set_xlim((1e-4,1e-2))
            #ax3.set_ylim([tmid.min(),tmid.max()])
            ax3.set_ylabel("$\\omega [rad.s^{-1}]$")
            plt.colorbar()
            #ax3.set_xticklabels([])
    
        return Pyy,frq,tmid
        
    def filt(self, cutoff_dt, btype='low', order=3, axis=-1):
        """
        Butterworth filter the time series
        
        Inputs:
            cutoff_dt - cuttoff period [seconds]
            btype - 'low' or 'high' or 'band'
        """
        
        if self.isequal==False and self.VERBOSE:
            print('Warning - time series is unequally spaced.\
                Use self.interp to interpolate onto an equal grid')
        
        if not btype == 'band':
            Wn = self.dt/cutoff_dt
        else:
            Wn = [self.dt/co for co in cutoff_dt]
            
        (b, a) = signal.butter(order, Wn, btype=btype, analog=0, output='ba')
        
        # filtfilt only likes to operate along the last axis
        ytmp = np.swapaxes(self.y,-1,axis)
        ytmp = signal.filtfilt(b, a, ytmp, axis=-1)

        return np.swapaxes(ytmp,axis,-1)
        #return signal.filtfilt(b, a, self.y, axis=axis)

    def filt_uneven(self, cutoff_dt, axis=-1):
        """
        Lowpass filter with gaps in the output
        
        Inputs:
            cutoff_dt - cuttoff period [seconds]
            btype - 'low' or 'high' or 'band'
        """
        

        ndim = self.y.ndim
        assert ndim <= 2
        
        if ndim == 1:
            return filt_gaps(self.y.data, self.y.mask, self.dt, cutoff_dt)
        else:
            # operate along the last axis
            ytmp = np.swapaxes(self.y,-1,axis)
            yf = np.zeros_like(ytmp)
            nk = ytmp.shape[0]
            for kk in range(nk):
                yf[kk,:] = filt_gaps(ytmp.data[kk,:], ytmp.mask[kk,:], self.dt, cutoff_dt)

            return np.swapaxes(yf,axis,-1)


    def godinfilt(self,filtwidths=[24,25]):
        """
        Apply the successive Godin-type filter to the data set
        
        The filter widths are defined by filtwidths (hours).
        
        For a 24-25-24 hour filter set filtwidths=[24,25] (default).
        """
        
        if self.isequal==False or self.dt != 3600.:
            # Puts the data onto an hourly matrix 
            print('Warning -- downsampling data to use Godin filter...')
            y, t = self._evenly_dist_data(self.y, self.tsec, 3600.)
        else:
            y, t = self.y, self.t
            
        ymean = self.running_mean(y, 3600., windowlength=filtwidths[0]*3600.)
        ymean = self.running_mean(ymean, 3600., windowlength=filtwidths[1]*3600.)
        ymean = self.running_mean(ymean, 3600., windowlength=filtwidths[0]*3600.)

        # Zero out masks
        ymean[ymean.mask] = 0

        # Ensure the output is on the same time grid
        tout, yout = self.copy_like(t, ymean).interp(self.t, method='cubic', axis=-1)
        yout = np.ma.MaskedArray(yout, self.y.mask)
        return self.copy_like(tout, yout)
        #self.y = ymean
            
           
    def interp(self,timein,method='linear',timeformat='%Y%m%d.%H%M%S',axis=-1):
        """
        Interpolate the data onto an equally spaced vector
        
        timein is either:
            (3x1 tuple) - (tstart,tend,dt)
                tstart and tend - string with format 'yyyymmdd.HHMMSS'
        or
            datetime vector
        
        method - method passed to interp1d
               - use 'nearest' to preserve masking in gap regions
        """
        
        # Create the time vector
        try:
            tstart=timein[0]
            tend=timein[1]
            dt=timein[2]
            tnew = othertime.TimeVector(tstart,tend,dt,timeformat=timeformat)
            tsec = self._get_tsec(tnew)
        except:
            tnew = self._set_time(timein)
            tsec = self._get_tsec(tnew)
            #dt = tsec[1] - tsec[0]
            dt,isequal = self._check_dt(tsec)
            #dt = (tnew[1]-tnew[0]).total_seconds()

        if method=='nearestmask':
            # Nearest neighbour doesn't use interp1d to preserve mask
            #self._evenly_dist_data(dt)
            #tnew, output = self.subset(tnew[0],tnew[-1])
            output = self._interp_nearest(tsec)

        else:

            #t = othertime.SecondsSince(tnew,basetime = self.basetime)
            # Don't include nan points
            if self.ndim > 1:
                # Interpolate multidimensional arrays without a mask
                F = interpolate.interp1d(self.tsec, self.y,\
		    kind=method,axis=axis,\
                    bounds_error=False,fill_value=0)

                output = F(tsec)

            else:

                #mask = np.isnan(self.y) == False
                if np.all(self.y.mask):
                    output = np.zeros(t.shape)
                    output[:] = np.nan
                else:

                    if np.all(self.y.mask == False) :
                        mask = np.isfinite(self.y).data
                    else:
                        mask = ~self.y.mask

                    F = interpolate.interp1d(\
		    	self.tsec[mask], self.y[mask],\
			kind=method, axis=axis,\
                        bounds_error=False, fill_value=0)

                    output = F(tsec)

        return tnew, output
        
        
    def tidefit(self,frqnames=None,basetime=None,axis=-1):
        """
        Perform a tidal harmonic fit to the data
        
        Returns the amp, phase, frequencies and fitted time series
        """
        
        # Get the tidal fruequencies
        if frqnames is None:
            # This returns the default frequencies from the uspectra class
            frq,frqnames = getTideFreq(Fin=None)
        else:
            frq,frqnames = getTideFreq(Fin=frqnames)
            
        # Call the uspectra method
        #U = uspectra(self.t, self.y, frq=frq, method='lsqfast',axis=axis)
        #amp,phs = U.phsamp(phsbase=basetime)
        #return amp, phs, frq, frqnames, U.invfft()

        amp, phs, mean = \
            harmonic_fit(self.t, self.y, frq,\
             mask=self.y.mask,\
             phsbase=basetime, axis=axis)

        yfit = harmonic_signal(self.t, amp, phs, mean, frq, phsbase=basetime, axis=axis)
        
        # compute an error
        yrms = rms(self.y - yfit.T)
        
        return amp, phs, frq, mean, yfit, yrms
        #return amp, phs, frq, frqnames, yfit
        
    def running_harmonic(self,omega,windowlength=3*86400.0,overlap=12*3600.0, plot=True):
        """
        Running harmonic fit of the time series at frequency, omega. 
        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        """
        
        # Make sure that omega is a list
        try:
            len(omega)
        except:
            omega=[omega]
        
        pt1,pt2 = window_index_time(self.t,windowlength,overlap)
        npt = len(pt1)
        tmid = []
        amp = np.zeros((npt,))
        phs = np.zeros((npt,))
        ymean = np.zeros((npt,))
        ii=-1
        for t1,t2 in zip(pt1,pt2):
            ii+=1
            # Perform the harmonic fit on the segment
            U = uspectra(self.t[t1:t2],self.y[t1:t2],frq=omega,method='lsqfast')
            # Reference the phase to the start of the time series
            amp[ii],phs[ii] = U.phsamp(phsbase = self.t[0])
            
            # Return the mid time point
            ind = np.floor(t1 + (t2-t1)/2)
            tmid.append(self.t[ind])
            
            # Return the fitted time series
            ymean[ii] = self.y[t1:t2].mean()
            #yout[t1:t2] += ymean + U.invfft()

            
        tmid = np.asarray(tmid)
        

        
        if plot:
            plt.subplot(211)
            self.plot()
            plt.plot(tmid,ymean,'r')
            plt.fill_between(tmid,ymean-amp,y2=ymean+amp,color=[0.5,0.5,0.5],alpha=0.5)
            plt.legend(('Original Signal','Harmonic reconstruction'))
            
            ax=plt.subplot(212)
            plt.fill_between(tmid,amp,alpha=0.5)
            plt.xticks(rotation=17)
            plt.ylabel('Amplitude')
            ax.set_xlim([self.t[0],self.t[-1]])
            
            
        return tmid, amp, phs
            
    def running_mean(self, y, dt, windowlength=3*86400.0):
        """
        Running mean of the time series
        
        windowlength - length of each time window [seconds]
        """
        mask = y.mask.copy()
        y[y.mask]=0.
        y.mask=mask
        
        windowsize = int(np.floor(windowlength/dt))
        ytmp = y.copy()
        ytmp = self._window_matrix(ytmp,windowsize)
        
        weights = 1./windowsize * np.ones((windowsize,))
        #ypt = ne.evaluate("ytmp * weights")
        #ypr = ytmp*weights[np.newaxis,...]
        #ypr = ytmp.__mul__(weights)
        #ytmp2 = ypr.sum(axis=-1)
        ytmp2 = np.sum(ytmp*weights,axis=-1)
        
        # This result needs to be normalized to account for missing data,
        # this is the same as calculating different weights for each section
        ntmp= np.ones_like(y)
        ntmp[mask] = 0.
        norm = self._window_matrix(ntmp,windowsize)
        #norm*= weights
        norm = norm.sum(axis=-1)
        norm /= windowsize
        
        ytmp2/=norm
                
        return self._update_windowed_data(y, ytmp2,windowsize)

    def running_rms(self, y, windowlength=3*86400.0):
        """
        Running RMS of the time series

        windowlength - length of each time window [seconds]
        """
        mask = y.mask.copy()
        y[y.mask]=0.
        y.mask=mask
        
        windowsize = int(np.floor(windowlength/self.dt))
        ytmp = y.copy()
        ytmp = self._window_matrix(ytmp,windowsize)
        ytmp2 = np.sum(ytmp*ytmp,axis=-1)

        ntmp= np.ones_like(self.y)
        ntmp[mask] = 0.
        N = self._window_matrix(ntmp,windowsize)
        N = N.sum(axis=-1)

        return self._update_windowed_data(y, np.sqrt( 1./N * ytmp2),windowsize)

    def despike(self,nstd=4.,windowlength=3*86400.0,overlap=12*3600.0,\
        upper=np.inf,lower=-np.inf,maxdiff=np.inf,fillval=0.):
        """
        Despike time series by replacing any values greater than nstd*std.dev with the
        median of a running window. 
        
        nstd - number of standard deviations outside to replace        
        windowlength - length of each time window [seconds]
        overlap - overlap between windows [seconds]
        lower - lower value bound
        upper - upper value bound
        maxdiff - maximum difference between two points
        """
#        if self.isequal==False:
#            self._evenly_dist_data()
            
        nbad = 0
        
        # Now check the maximum difference
        ydiff = np.zeros_like(self.y)
        ydiff[1::] = np.abs(self.y[1::]-self.y[0:-1])
        ind = ydiff>=maxdiff        
        #self.y[ind]=fillval
        self.y.mask[ind] = True # Mask needs to be set after values are prescribed
        nbad += np.sum(ind)
        
        # First mask any NaN and values outside of bounds
        ind = operator.or_(self.y<=lower,self.y>=upper)
        #self.y[ind]=fillval
        self.y.mask[ind] = True
        nbad += np.sum(ind)
        
        ind =  np.isnan(self.y)
        #self.y[ind]=fillval
        self.y.mask[ind] = True
        
        nbad += np.sum(ind)
        
        # Now calculate the moving median and standard deviation
        windowsize = np.floor(windowlength/self.dt)
        ytmp = self.y.copy()
        ytmp = self._window_matrix(ytmp,windowsize)
        
        ytmp2 = np.mean(ytmp,axis=-1)
        ymean = self._update_windowed_data(self.y, ytmp2,windowsize)
        
        #ytmp2= np.std(ytmp,axis=-1)
        ytmp2 = np.apply_along_axis(np.std,-1,ytmp2)
        ystd = self._update_windowed_data(self.y, ytmp2,windowsize)
        
        # Mask values outsize of the
        ind = operator.or_(self.y >= ymean + nstd*ystd,\
                self.y <= ymean - nstd*ystd)
        
        #self.y[ind] = ymedian[ind]
        self.y.mask[ind] = True
        
        nbad += np.sum(ind)
        
        if self.VERBOSE:
            print('Despiked %d points'%nbad)
        
        
        
    def find_trend(self, axis=-1):
        """
        Return the linear trend as a numpy array
        """
        # Make sure time is along the last axis
        ytmp = np.swapaxes(self.y,-1,axis)

        # Make sure there are no NaNs
        ytmp[ytmp.mask] = 0.0

        # Fit
        model = np.polyfit(self.tsec, ytmp.data.T, 1) # Linear fit

        # 
        Nz = model.shape[1]
        trend = np.zeros((Nz, self.Nt))
        for kk in range(Nz):
             trend[kk,:] = np.polyval(model[:,kk], self.tsec)

        return trend

    def plot(self,angle=17,**kwargs):
        """
        Plot
        
        Rotates date labels
        """
        
        h1=plt.plot(self.t,self.y.T,**kwargs)
        plt.xticks(rotation=angle)
        
        return h1 
        
    def fillplot(self,angle=17,alpha=0.7,**kwargs):
        """
        
        """
        h1=plt.fill_between(self.t,self.y,alpha=alpha,**kwargs)
        plt.xticks(rotation=angle)
        
        return h1 

    def get_tslice(self,time1,time2):
        """
        Returns the time indices bounded by time1 and time2
        """
        t0 = othertime.findNearest(time1,self.t)
        t1 = othertime.findNearest(time2,self.t)
        #try:
        #    t0 = othertime.findNearest(time1,self.t)
        #    t1 = othertime.findNearest(time2,self.t)
        #except:
        #    t0 = np.searchsorted(self.t, np.datetime64(time1))
        #    t1 = np.searchsorted(self.t, np.datetime64(time2))

        #t1 = min( self.Nt, t1+1)
        t1+=1
        if t1 > self.t.shape[0]:
            t1=self.t.shape[0]
        #t0-=1
        #if t0<0:
        #    t0 = 0

        return t0, t1

    def clip(self, t0, t1):
        """
        Clip the object between the time limits t1 - t2

        Differs from subset as it returns a timeseries object

        Inputs:
                t1, t2 : datetime
        Returns:
                a NEW timeseries object
        """
        #t0 = datetime.strptime(t0, fmt)
        #t1 = datetime.strptime(t1, fmt)

        modt,mody = self.subset(t0,t1)

        return self.copy_like(modt, mody)

        
    def subset(self,time1,time2):
        """
        Returns a subset of the array between time1 and time2
        """
        t0, t1 = self.get_tslice(time1, time2)
           
        return self.t[t0:t1], self.y[..., t0:t1]
    
    def resample(self, dt, ndt=2.):
        """
        Resamples the current timeseries object at the output interval
        "dt"

        Applies a running mean (window=dt) if dt > ndt*self.dt
        """

        # Output time vector
        timeout = othertime.TimeVector(self.t[0], self.t[-1], dt, 
            istimestr=False)

        # Use a running mean to filter data
        if dt > ndt*self.dt:
            ymean = self.running_mean(self.y, self.dt, windowlength=dt)
            # Create a time series with the filtered data
            tsmean = self.copy_like(self.t, ymean)

            # Interpolate onto the output step
            tnew, ynew = tsmean.interp(timeout, method='nearestmask')

            return self.copy_like(tnew, ynew)
        
        else:
            # Do not filter if dt < self.dt
            # Interpolate onto the output step
            tnew, ynew = self.interp(timeout, method='linear')

            return self.copy_like(tnew, ynew)

	
	        
    def copy_like(self, t, y):
        """
        Create a new time series "like" current one

        Retains relevant attributes
        """

        return timeseries(t, y,\
            units='%s'%self.units,\
            long_name='%s'%self.long_name,\
            StationID = '%s'%self.StationID,\
            StationName = '%s'%self.StationName,\
            varname = '%s'%self.varname,\
            X= self.X,\
            Y= self.Y,\
            Z= self.Z,\
        )

    def savetxt(self,txtfile):
        f = open(txtfile,'w')
        
        #for t,v in zip(self.tsec,self.y):
        #    f.write('%10.6f\t%10.6f\n'%(t,v))
        for ii in range(self.y.shape[0]):
            f.write('%s, %10.6f\n'%(datetime.strftime(self.t[ii],'%Y-%m-%d %H:%M:%S'),self.y[ii]))
            
        f.close()

    def to_xray(self):
        """
        Convert the timeseries to an xray.DataArray object
        """
        # Define the coordinates
        coords = {'time':self.t.copy()}
        if self.y.ndim == 1:
            dims = ('time',)
            y = self.y
        else:
            dims = ('time','depth')
            #dims = ('depth','time')
            coords.update({'depth':self.Z})
            y = self.y.T 

        # Define the attributes
        attrs = {\
            'units':'%s'%self.units,\
            'long_name':'%s'%self.long_name,\
            'StationID':'%s'%self.StationID,\
            'StationName':'%s'%self.StationName,\
        }

        return xray.DataArray( y.copy(), \
                    dims=dims,\
                    name=self.varname,\
                    attrs = attrs,\
                    coords = coords,\
               )
 
    def copy(self):
        """
        Make a copy of the time-series object in memory
        """
        from copy import deepcopy
        return deepcopy(self)
        
    def _get_tsec(self, time):

        # Convert the time to seconds
        #if isinstance(time[0], np.datetime64):
        #    tsec = ((time - time[0])*1e-9).astype(np.float64)
        #else:
        #    tsec = othertime.SecondsSince(time,basetime=self.basetime)

        # SecondsSince handle datetime64 objects
        tsec = othertime.SecondsSince(time,basetime=self.basetime)

        return tsec
 
    def _set_time(self, t):
        """
        Return the time variable as a datetime64 object
        ~~Return the time variable as a datetime object~~
        """
        if isinstance(t, pd.DatetimeIndex):
            #return othertime.datetime64todatetime(t.values)
                return t.values
        elif isinstance(t[0], datetime):
            return othertime.datetimetodatetime64(t) 
        elif isinstance(t[0], np.datetime64):
                return t
            #return othertime.datetime64todatetime(t) 
        else:
            raise Exception('unknown time type: ').with_traceback(type(t))

        #if isinstance(t, pd.tseries.index.DatetimeIndex):
        #    return othertime.datetime64todatetime(t.values)
        #elif isinstance(t[0], datetime):
        #    return t
        #elif isinstance(t[0], np.datetime64):
        #    return othertime.datetime64todatetime(t) 
        #else:
        #    raise Exception, 'unknown time type: ', type(t)

    def _check_dt(self, tsec):
        """
        Check that the time series is equally spaced
        """
        dt = np.diff(tsec)
        
        dt_unique = np.unique(dt)
        
        if np.size(dt_unique) == 1:
            isequal = True
        else:
            isequal = False
        
        try:
            dt = dt[1]
        except:
            dt = 0.0

        return dt, isequal
            
    def _interp_nearest(self, tout, y=None):
        """
        Nearest interpolation with preseverd mask
        """
        if y is None:
            y = self.y

        t0 = self.tsec[0]
        t = self.tsec - t0

        #tout -= tout[0]
        tnew = tout - t0 # Reference to the same time

        shape = y.shape[:-1] + tnew.shape
        yout = np.ma.MaskedArray(np.zeros(shape), mask=True)

        tind = np.searchsorted(tnew, t) - 1

        yout[...,tind] = y[:]

        return yout

    def _evenly_dist_data(self, y, tsec, dt):
        """
        Distribute the data onto an evenly spaced array
        
        No interpolation is performed
        """
        if self.VERBOSE:
            print('inserting the data into an equally-spaced time vector (dt = %f s).'%self.dt)
    
        t0 = tsec[0]
        t = tsec - t0
        # Put the data onto an evenly spaced, masked array
        tout = np.arange(t[0],t[-1]+dt,dt)
        
        tind = np.searchsorted(tout,t) - 1
        
        shape = y.shape[:-1] + tout.shape
        yout = np.ma.MaskedArray(np.zeros(shape),mask=True)

        yout[...,tind] = y

        ### Slow method
        #def updatetime(tsec):
        #    try:
        #        return np.timedelta64(int(tsec), 's') + self.t[0]
        #    except:
        #        return timedelta(seconds=tsec) + self.t[0]
        #    
        #self.t = np.array(map(updatetime,tout))

        t = tout.astype("timedelta64[s]") + np.datetime64(self.t[0])

        #self.y = yout
        
        #self.tsec = tout+t0
        #
        #self.ny = np.size(y)
        #
        #self.isequal = True
        #
        #self.dt = dt

        return yout, t
        
    def _window_matrix(self,y,windowsize):
        """
        Returns the matrix as a strided array so that 'rolling' operations can
        be performed along the last axis
        """
        windowsize=int(windowsize)
        shape = y.shape[:-1] + (y.shape[-1] - windowsize + 1, windowsize)
        strides = y.strides + (y.strides[-1],)
        
        # The masked values get 
        return np.lib.stride_tricks.as_strided(y, shape=shape, strides=strides)
    
    def _update_windowed_data(self, y, ytmp,windowsize):
        """
        Re-inserts data that has been windowed into an array
        that is the same size as the original time series
        """
        y0 = np.zeros_like(y)
        indent = (windowsize-np.mod(windowsize,2))/2
        
        if np.mod(windowsize,2)==1:
            y0[...,indent:-indent]=ytmp
        else:
            y0[...,indent-1:-indent]=ytmp
        
        y0 = np.ma.MaskedArray(y0,mask=y.mask)
        y0.mask[...,0:indent]=True
        y0.mask[...,-indent:]=True
        
        return y0

class ModVsObs(object):
    """
    Class for handling and comparing two time series i.e. model vs observations

    DEPRECATED - moved to its own module modvsobs.py
    """

    units=' '
    long_name=' '
    StationID = ' '
    varname = ' '
    Z=0.0

    def __init__(self,tmod,ymod,tobs,yobs,**kwargs):
        """
        Inputs:
            tmod,tobs - vector of datetime object
            ymod,yobs - vector of values

        Keywords:
            long_name: string containing variable's name (used for plotting)
            units: string containing variable's units (used for plotting)

        Note that tmod and tobs don't need to be the same length. yobs is
        linearly interpolated onto tmod.
        """
        self.__dict__.update(kwargs)


        # Set the range inclusive of both observation and model result
        time0 = max(tmod[0],tobs[0])
        time1 = min(tmod[-1],tobs[-1])

        if time1 < time0:
            print('Error - the two datasets have no overlapping period.')
            return None
        
        # Clip both the model and observation to this daterange

        t0 = othertime.findNearest(time0,tmod)
        t1 = othertime.findNearest(time1,tmod)
        TSmod = timeseries(tmod[t0:t1],ymod[...,t0:t1])

        t0 = othertime.findNearest(time0,tobs)
        t1 = othertime.findNearest(time1,tobs)
        self.TSobs = timeseries(tobs[t0:t1],yobs[...,t0:t1])

        # Interpolate the observed value onto the model step
        #tobs_i, yobs_i = TSobs.interp(tmod[t0:t1],axis=0)
        #self.TSobs = timeseries(tobs_i, yobs_i)

        # Interpolate the modeled value onto the observation time step
        tmod_i, ymod_i = TSmod.interp(tobs[t0:t1],axis=-1)
        self.TSmod = timeseries(tmod_i,ymod_i)

        self.N = self.TSmod.t.shape[0]
        if self.N==0:
            print('Error - zero model points detected')
            return None

        self.calcStats()

        # Calculate the data limits
        self.calc_data_lims()

    def calc_data_lims(self):
        y0 = min(self.TSobs.y.min(), self.TSmod.y.min())
        y1 = max(self.TSobs.y.max(), self.TSmod.y.max())
        #ymax = max[np.abs(y0), np.abs(y1)]
        self.ylims = [y0, y1]

    def plot(self,colormod='r',colorobs='b',legend=True,loc='lower right',**kwargs):
        """
        Time-series plots of both data sets with labels
        """
        ax = plt.gca()

        h1 = self.TSmod.plot(color=colormod,**kwargs)

        h2 = plt.plot(self.TSobs.t,self.TSobs.y,color=colorobs,**kwargs)

        plt.ylabel(r'%s [$%s$]'%(self.long_name,self.units)) # Latex formatting

        plt.grid(b=True)

        plt.title('StationID: %s'%self.StationID)

        if legend:
            plt.legend(('Model','Observed'),loc=loc)

        return h1, h2, ax
        
    def stackplot(self,colormod='r',colorobs='b',scale=None,ax=None,fig=None,labels=True,**kwargs):
        """
        Stack plot of several time series
        """
        if labels:
            labels = ['z = %1.1f m'%z for z in self.Z.tolist()]
        else:
            labels=None
        
        fig,ax,ll = stackplot(self.TSobs.t,self.TSobs.y,ax=ax,fig=fig,\
            scale=scale,units=self.units,labels=labels,color=colorobs,*kwargs)
            
        fig,ax,ll = stackplot(self.TSmod.t,self.TSmod.y,ax=ax,fig=fig,\
            scale=scale,units=self.units,labels=labels,color=colormod,*kwargs)

    def scatter(self,ylims=None,printstats=True,\
           textw=0.05, texth=0.55, **kwargs):
        """
        Scatter plot of the model vs observation
        """
        if ylims==None:
            ylims = self.ylims

        h1 = plt.plot(self.TSobs.y.ravel(),self.TSmod.y.ravel(),'.',**kwargs)
        plt.plot([ylims[0],ylims[1]],[ylims[0],ylims[1]],'k--')

        ax = plt.gca()
        ax.set_aspect('equal')
        plt.xlim(ylims)
        plt.ylim(ylims)
        #ax.autoscale(tight=True)
        plt.grid(b=True)

        if printstats:
            textstr = '$r^2$ = %6.2f\nRMSE = %6.2f\nBias = %6.2f\n'%(\
                self.cc.mean(),self.rmse.mean(),self.bias.mean())
            plt.text(textw,texth,textstr,transform=ax.transAxes)

        return h1, ax

    def qqplot(self, percentiles=[1.,5.,25.,50.,75.,95.,99.],\
                ylims=None, **kwargs):
        """
        Quantile-quantile plot
        """
        q_mod = np.percentile(self.TSmod.y, percentiles)
        q_obs = np.percentile(self.TSobs.y, percentiles)

        if ylims==None:
            ylims = self.ylims

        # scale the marker size
        sizes = (1 - np.abs(np.array(percentiles)-50)/50)*50


        h1 = plt.scatter(q_obs, q_mod,s=sizes, **kwargs)
        plt.plot([ylims[0],ylims[1]],[ylims[0],ylims[1]],'k--')

        ax = plt.gca()
        ax.set_aspect('equal')
        plt.xlim(ylims)
        plt.ylim(ylims)
        #ax.autoscale(tight=True)
        plt.grid(b=True)

        return h1, ax

    def calcStats(self):
        """
        Calculates statistics including:
            moments, RMS, CC, skill, ...
        """
        self.meanObs = self.TSobs.y.mean(axis=-1)
        self.meanMod = self.TSmod.y.mean(axis=-1)
        self.stdObs = self.TSobs.y.std(axis=-1)
        self.stdMod = self.TSmod.y.std(axis=-1)
        self.rmsObs = rms(self.TSobs.y,axis=-1)
        self.rmsMod = rms(self.TSmod.y,axis=-1)

        # RMSE
        self.rmse = rms(self.TSobs.y-self.TSmod.y,axis=-1)
        
        # bian
        self.bias = np.mean(self.TSmod.y - self.TSobs.y, axis=-1)

        # skill
        self.skill = 1.0 - ((self.TSobs.y-self.TSmod.y)**2.).sum(axis=-1) / \
            ( (self.TSobs.y.T - self.meanObs)**2.).T.sum(axis=-1)

        # Correlation coefficient
        self.cc = 1.0/float(self.N) * ( (self.TSobs.y.T-self.meanObs).T * \
            (self.TSmod.y.T - self.meanMod).T ).sum(axis=-1) / (self.stdObs * self.stdMod) 

    def printStats(self,f=None,header=True):
        """
        Prints the statistics to a markdown language style table
        """
        if 'meanMod' not in self.__dict__:
            self.calcStats()

        outstr=''

        if header:
            outstr += "|      | Mean Model | Mean Obs. | Std. Mod. | Std Obs | RMSE |   CC   | skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        outstr += "| %s [%s] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (self.StationID,self.units, self.meanMod, self.meanObs, self.stdMod, self.stdObs,\
            self.rmse,self.cc,self.skill)

        if f == None:
            print(outstr)
        else:
            f.write(outstr)


    def printStatsZ(self,f=None,header=True):
        """
        Prints the statistics to a markdown language style table
        """
        outstr=''

        if header:
            outstr += "| Depth | Mean Model | Mean Obs. | Std. Mod. | Std Obs | RMSE |   CC   | skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        for ii,zz in enumerate(self.Z.tolist()):

            outstr += "| %3.1f [m] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (zz, self.meanMod[ii], self.meanObs[ii], self.stdMod[ii],\
                self.stdObs[ii], self.rmse[ii],self.cc[ii],self.skill[ii])

        if f == None:
            print(outstr)
        else:
            f.write(outstr)

    def crosscorr(self,normalize=False,axis=-1):
        """
        Crosscorrelation calculation
        """

        assert self.TSobs.isequal,\
             'Data must be equally spaced to perform this function.'

        N = self.TSobs.ny
        M = int(N)/10

        ymean = self.TSobs.y.mean(axis=axis)
        y = self.TSobs.y - ymean
        xmean = self.TSmod.y.mean(axis=axis)
        x = self.TSmod.y - xmean

        k = list(range(1,M))
        tau = np.asarray(k,dtype=np.float)*self.TSobs.dt

        Cxy = [1./(N-kk) * np.sum(y[...,0:-kk]*x[...,kk::],axis=axis) for kk in k ]

        if normalize:
            return Cxy/(y.std()*x.std()), tau
        else:
            return Cxy ,tau
 
    def csd(self, plot=True,nbandavg=1,**kwargs):
        """
        Cross spectral density
        
        nbandavg = Number of fft bins to average
        """
        
        if self.isequal==False and self.VERBOSE:
            print('Warning - time series is unequally spaced consider using lomb-scargle fft')
        

        NFFT = int(2**(np.floor(np.log2(self.ny/nbandavg)))) # Nearest power of 2 to length of data
            
        Pyy,frq = mlab.csd(self.TSobs.y-self.TSobs.y.mean(),Fs=2*np.pi/self.dt,NFFT=NFFT,window=mlab.window_hanning,scale_by_freq=True)
        
        if plot:
            plt.loglog(frq,Pyy,**kwargs)
            plt.xlabel('Freq. [$cycles s^{-1}$]')
            plt.ylabel('PSD')
            plt.grid(b=1)
        
        return Pyy, frq


def loadDBstation(dbfile, stationName, varname, timeinfo=None, \
     filttype=None,cutoff=3600.0,output_meta=False,method='linear'):
    """
    Load station data from a database file
    
    Inputs:
        dbfile - location of database file
        stationName - StationName in database
        varname - variable name e.g. 'waterlevel', 'discharge', 'salinity'
        
        timeinfo (optional) - tuple with (starttime,endtime,dt). Format 'yyyymmdd.HHMMSS'
            Use this to interpolate onto a constant time vector
        filttype (optional) - 'low' or 'high' 
            Set this to filter data
            
    Returns:
        timeseries object
        -1 on error
            
    """
    
    outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName']
    tablename = 'observations'
    condition = 'Variable_Name = "%s" and StationID = "%s"' % (varname,stationName)
    #condition = 'Variable_Name = "%s" and StationName LIKE "%%%s%%"' % (varname,stationName)
    
    print('Querying database...')
    print(condition)
    data, query = queryNC(dbfile,outvar,tablename,condition)  

    yout = data[0][varname].squeeze()
    # Zero nan
    yout[np.isnan(yout)] = 0.0
    
    if len(data)==0:
        print('!!! Warning - Did not find any stations matching query. Returning -1 !!!')
        return -1
    else:
        ts = timeseries(data[0]['time'],yout)
        
        
    if not timeinfo==None:
        print('Interpolating station data between %s and %s\n'%(timeinfo[0],timeinfo[1]))
        tnew,ynew =\
            ts.interp((timeinfo[0],timeinfo[1],timeinfo[2]),method=method)
        ts = timeseries(tnew,ynew)
        ts.dt = timeinfo[2] # This needs updating
        
    if not filttype==None:
        print('%s-pass filtering output data. Cutoff period = %f [s].'%(filttype,cutoff))
        yfilt = ts.filt(cutoff,btype=filttype,axis=-1)
        ts.y = yfilt.copy()
    
    if output_meta:
        if 'elevation' in data[0]:
            ele = data[0]['elevation']
        else:
            ele = np.array([0.0])
        meta = {'longitude':data[0]['longitude'],'latitude':data[0]['latitude'],'elevation':ele,'StationName':query['StationName'][0]}
        return ts, meta        
    else:
        return ts

def window_index_time_slow(t,windowsize,overlap):
    """
    Determines the indices for window start and end points of a time vector
    
    The window does not need to be evenly spaced
    
    Inputs:
        t - list or array of datetime objects
        windowsize - length of the window [seconds]
        overlap - number of overlap points [seconds]
        
    Returns: pt1,pt2 the start and end indices of each window
    """
    
    try:
        t=t.tolist()
    except:
        t=t
        
    t1=t[0]
    t2=t1 + timedelta(seconds=windowsize)
    pt1=[0]
    pt2=[othertime.findNearest(t2,t)]
    while t2 < t[-1]:
        t1 = t2 - timedelta(seconds=overlap)
        t2 = t1 + timedelta(seconds=windowsize)

        pt1.append(othertime.findNearest(t1,t))
        pt2.append(othertime.findNearest(t2,t))
        
    return pt1, pt2
    
def window_index_time(t,windowsize,overlap):
    """
    Determines the indices for window start and end points of a time vector
    
    The window does not need to be evenly spaced
    
    Inputs:
        t - list or array of datetime objects
        windowsize - length of the window [seconds]
        overlap - number of overlap points [seconds]
        
    Returns: pt1,pt2 the start and end indices of each window
    """
    
    tsec = othertime.SecondsSince(t)
        
    t1=tsec[0]
    t2=t1 + windowsize
    pt1=[0]
    pt2=[np.searchsorted(tsec,t2)]
    while t2 < tsec[-1]:
        t1 = t2 - overlap
        t2 = t1 + windowsize

        pt1.append(np.searchsorted(tsec,t1))
        pt2.append(np.searchsorted(tsec,t2))
        
    return pt1, pt2
    
def pol2cart(th,rho):
    """Convert polar coordinates to cartesian"""
    x = rho * np.cos(th)
    y = rho * np.sin(th)

    return x, y
    
def cart2pol(x,y):
    """
    Convert cartesian to polar coordinates
    """
    th = np.angle(x+1j*y)
    rho = np.abs(x+1j*y)
    
    return th, rho
    
def ap2ep(uamp,uphs,vamp,vphs):
    """
    Convert u/v amplitude phase information to tidal ellipses
    
    All angles are in radians
    
    Returns:
        SEMA, SEMI, INC, PHS, ECC
    
    Based on the MATLAB ap2ep function:
        https://www.mathworks.com/matlabcentral/fileexchange/347-tidalellipse/content/ap2ep.m
    """
    # Make complex amplitudes for u and v
    u = uamp*np.exp(-1j*uphs)
    v = vamp*np.exp(-1j*vphs)

    #Calculate complex radius of anticlockwise and clockwise circles:
    wp = (u+1j*v)/2.0     # for anticlockwise circles
    wm = np.conj(u-1j*v)/2.0  # for clockwise circles
    # and their amplitudes and angles
    Wp = np.abs(wp)
    Wm = np.abs(wm)
    THETAp = np.angle(wp)
    THETAm = np.angle(wm)
   
    # calculate ep-parameters (ellipse parameters)
    SEMA = Wp+Wm              # Semi  Major Axis, or maximum speed
    SEMI = Wp-Wm            # Semin Minor Axis, or minimum speed
    ECC = SEMI/SEMA          # Eccentricity

    PHA = (THETAm-THETAp)/2.0   # Phase angle, the time (in angle) when 
                               # the velocity reaches the maximum
    INC = (THETAm+THETAp)/2.0   # Inclination, the angle between the 
                               # semi major axis and x-axis (or u-axis).
                               
    return SEMA, SEMI, INC, PHA, ECC
    
def rms(x,axis=-1):
    """
    root mean squared
    """
    
    return np.sqrt(1.0/np.shape(x)[-1] * np.sum(x**2,axis=axis))
    
def crms(t,y):
    """
    Continous function rms
    """
    fac = 1.0/(t[-1]-t[0])
    return np.sqrt(fac*np.trapz(y**2,x=t))

def rmse(x1, x2, axis=-1):
    """
    Root mean squared error
    """
    return rms(x1-x2, axis=axis)

def skill(xmod, xobs, axis=-1):
    """
    Murphy skill score
    """
    varmod = ((xobs - xmod)**2.).sum(axis=axis)
    varobs = ((xobs - xobs.mean(axis=axis))**2.).sum(axis=axis)

    return 1 - varmod/varobs
    
def tidalrmse(Ao,Am,Go,Gm):
    """
    Tidal harmonic RMSE
    Ao, Am - observed and modeled amplitude
    Go, Gm - observed and modeled phase (radians)
    """
    return np.sqrt( 0.5*(Ao**2 + Am**2) - Ao*Am*np.cos(Go-Gm) )

def loadtxt(txtfile):
    """
    Loads a text file with two columns
    Column 1 is time: seconds since 1990-1-1
    Column 2 is the data.
    """
    f = open(txtfile,'r')
    
    t=[]
    y=[]
    for line in f.readlines():
        line = line.strip()
        ll = line.split(',')
        t.append(datetime.strptime(ll[0],'%Y-%m-%d %H:%M:%S'))
        y.append(float(ll[1]))
        
    f.close()
    return timeseries(np.array(t),np.array(y))


