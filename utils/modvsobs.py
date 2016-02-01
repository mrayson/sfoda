"""
Model vs Observations modules
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from matplotlib.lines import Line2D
from scipy import signal, interpolate
import othertime
from datetime import datetime, timedelta
import xray

from soda.utils.timeseries import timeseries, rms

import pdb

class ModVsObs(object):
    """
    Class for handling and comparing two time series i.e. model vs observations

    """

    units=' '
    long_name=' '
    stationid = ' '
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
        self.__dict__.update(**kwargs)


        # Set the range inclusive of both observation and model result
        time0 = max(tmod[0],tobs[0])
        time1 = min(tmod[-1],tobs[-1])

        if time1 < time0:
            print 'Error - the two datasets have no overlapping period.'
            return None
        
        # Clip both the model and observation to this daterange

        t0 = othertime.findNearest(time0,tmod)
        t1 = othertime.findNearest(time1,tmod)
        TSmod = timeseries(tmod[t0:t1],ymod[...,t0:t1], **kwargs)

        t0 = othertime.findNearest(time0,tobs)
        t1 = othertime.findNearest(time1,tobs)
        self.TSobs = timeseries(tobs[t0:t1],yobs[...,t0:t1], **kwargs)

        # Interpolate the observed value onto the model step
        #tobs_i, yobs_i = TSobs.interp(tmod[t0:t1],axis=0)
        #self.TSobs = timeseries(tobs_i, yobs_i)

        # Interpolate the modeled value onto the observation time step
        tmod_i, ymod_i = TSmod.interp(tobs[t0:t1],axis=-1)
        self.TSmod = timeseries(tmod_i,ymod_i)

        self.N = self.TSmod.t.shape[0]
        if self.N==0:
            print 'Error - zero model points detected'
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

        plt.title('StationID: %s'%self.stationid)

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
        if not self.__dict__.has_key('meanMod'):
            self.calcStats()

        outstr=''

        if header:
            outstr += "|      | Mean Model | Mean Obs. | Std. Mod. | Std Obs | RMSE |   CC   | skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        outstr += "| %s [%s] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (self.stationid,self.units, self.meanMod, self.meanObs, self.stdMod, self.stdObs,\
            self.rmse,self.cc,self.skill)

        if f == None:
            print outstr
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
            print outstr
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

        k = range(1,M)
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
            print 'Warning - time series is unequally spaced consider using lomb-scargle fft'
        

        NFFT = int(2**(np.floor(np.log2(self.ny/nbandavg)))) # Nearest power of 2 to length of data
            
        Pyy,frq = mlab.csd(self.TSobs.y-self.TSobs.y.mean(),Fs=2*np.pi/self.dt,NFFT=NFFT,window=mlab.window_hanning,scale_by_freq=True)
        
        if plot:
            plt.loglog(frq,Pyy,**kwargs)
            plt.xlabel('Freq. [$cycles s^{-1}$]')
            plt.ylabel('PSD')
            plt.grid(b=1)
        
        return Pyy, frq

    ###############
    # IO routines
    ###############
    def to_xray(self, attrs={}):
        """
        Converts both time series to an xray dataset object

        attrs = global attribute dictionary
        """

        OBS = self.TSobs.to_xray()
        MOD = self.TSmod.to_xray()

        varobs = '%s_obs'%self.varname
        varmod = '%s_mod'%self.varname

        return xray.Dataset({
                varobs:OBS,\
                varmod:MOD,\
            }, attrs=attrs)
                
    def to_netcdf(self, ncfile, ncgroup=None, mode='w', attrs={}):
        """
        Save the ModVsObs object to a netcdf-4 file

        Overwrite any old files by detault
        """

        ds = self.to_xray(attrs=attrs)

        ds.to_netcdf(ncfile, group=ncgroup, format='NETCDF4', mode=mode)
