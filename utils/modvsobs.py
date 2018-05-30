"""
Model vs Observations modules

Use function:
    load_netcdf(ncfile, group)

To load from a netcdf file
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import mlab
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from scipy import signal, interpolate
from datetime import datetime, timedelta
import xray
import operator


from soda.utils import othertime
from soda.utils.otherplot import stackplot
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

    def __init__(self,tmod,ymod,tobs,yobs, interpmodel=True, **kwargs):
        """
        Inputs:
            tmod,tobs - vector of datetime object
            ymod,yobs - vector of values
            interpmodel - [default: True] interp the model onto the observations

        Keywords:
            long_name: string containing variable's name (used for plotting)
            units: string containing variable's units (used for plotting)

        Note that tmod and tobs don't need to be the same length. yobs is
        linearly interpolated onto tmod.
        """
        self.__dict__.update(**kwargs)


        # Set the range inclusive of both observation and model result
        #if isinstance(tmod,list):
        #    time0 = max(tmod[0],tobs[0])
        #    time1 = min(tmod[-1],tobs[-1])
        #elif isinstance(tmod[0], np.datetime64):
        #    time0 = max(tmod[0],tobs[0])
        #    time1 = min(tmod[-1],tobs[-1])
	time0 = max(tmod[0],tobs[0])
	time1 = min(tmod[-1],tobs[-1])

        if time1 < time0:
            print('Error - the two datasets have no overlapping period.')
            return None
        
	if not (tmod.shape[0] == tobs.shape[0]) and\
		not (tmod[0] == tobs[0]) and not (tmod[-1] == tobs[-1]) :
	    # Clip both the model and observation to this daterange
	    t0m = othertime.findNearest(time0,tmod)
	    t1m = othertime.findNearest(time1,tmod)
	    TSmod = timeseries(tmod[t0m:t1m],ymod[...,t0m:t1m], **kwargs)

	    t0 = othertime.findNearest(time0,tobs)
	    t1 = othertime.findNearest(time1,tobs)
	    TSobs = timeseries(tobs[t0:t1],yobs[...,t0:t1], **kwargs)

	    # Interpolate the observed value onto the model step
	    #tobs_i, yobs_i = TSobs.interp(tmod[t0:t1],axis=0)
	    #self.TSobs = timeseries(tobs_i, yobs_i)

	    ## Don't interpolate if datasets are the same
	    #if np.all(tobs==tmod):
	    #   self.TSobs = TSobs
	    #   self.TSmod = TSmod
	    # Interpolate the modeled value onto the observation time step
	    if interpmodel:
		tmod_i, ymod_i = TSmod.interp(tobs[t0:t1],axis=-1,method='nearestmask')
		#self.TSmod = timeseries(tmod_i,ymod_i, **kwargs)
		self.TSmod = timeseries(tobs[t0:t1], ymod_i, **kwargs)
		self.TSobs = TSobs
	    else:
		tobs_i, yobs_i = TSobs.interp(tmod[t0m:t1m],axis=-1,method='nearestmask')
		#self.TSobs = timeseries(tobs_i,yobs_i, **kwargs)
		self.TSobs = timeseries(tmod[t0m:t1m], yobs_i, **kwargs)
		self.TSmod = TSmod
	else:

	    self.TSmod = timeseries(tmod, ymod, **kwargs)
	    self.TSobs = timeseries(tobs, yobs, **kwargs)



        ### Check the dimension sizes
        #print self.TSmod.y.shape, self.TSobs.y.shape
        #print self.TSmod.t.shape[0], self.TSobs.t.shape[0]
        assert self.TSmod.t.shape[0] == self.TSobs.t.shape[0],\
                'Number of time records not equal'

        assert self.TSmod.y.shape == self.TSobs.y.shape,\
                'Dimensions sizes not equal'



        self.N = self.TSmod.t.shape[0]
        if self.N==0:
            print('Error - zero model points detected')
            return None

        # Compute the error 
        self.error = self.TSmod.y - self.TSobs.y

        self.calcStats()

        # Calculate the data limits
        self._calc_data_lims()

    def _calc_data_lims(self):
        y0 = min(self.TSobs.y.min(), self.TSmod.y.min())
        y1 = max(self.TSobs.y.max(), self.TSmod.y.max())
        #ymax = max[np.abs(y0), np.abs(y1)]
        self.ylims = [y0, y1]

    def clip(self, t0, t1, fmt='%Y-%m-%d'):
        """
        Clip the object between the time limits t1 - t2

        Inputs:
                t1, t2 : date string
                fmt : time format [default: '%Y-%m-%d'] 
        Returns:
                a NEW ModVsObs object
        """
        t0 = datetime.strptime(t0, fmt)
        t1 = datetime.strptime(t1, fmt)

        modt,mody = self.TSmod.subset(t0,t1)
        obst,obsy = self.TSobs.subset(t0,t1)

        return ModVsObs(modt, mody, obst, obsy,\
            units=self.units,\
            long_name=self.long_name,\
            stationid = self.stationid,\
            varname = self.varname,\
            Z= self.Z,\
        )



    def plot(self, colormod='r', colorobs='b', ylims=None, \
            legend=True, loc='lower right',
            dateformat=None,\
            **kwargs):
        """
        Time-series plots of both data sets with labels
        """
        ax = plt.gca()

        h1 = self.TSmod.plot(color=colormod,**kwargs)

        h2 = plt.plot(self.TSobs.t, self.TSobs.y.T, color=colorobs,**kwargs)

        if ylims is None:
            ylims = self.ylims

        ax.set_ylim(ylims)

        if dateformat is not None:
            ax.format_xdata = mdates.DateFormatter(dateformat)

        plt.ylabel(r'%s [$%s$]'%(self.long_name,self.units)) # Latex formatting

        plt.grid(b=True)

        plt.title('StationID: %s'%self.stationid)

        if legend:
            leg = plt.legend(('Model','Observed'),loc=loc)
            leg.get_frame().set_alpha(0.5)

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
        if ylims is None:
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
            textstr = 'skill = %6.2f\n$r^2$ = %6.2f\nRMSE = %6.2f\nBias = %6.2f\n'%(\
                self.skill.mean(),self.cc.mean(),self.rmse.mean(),self.bias.mean())
            plt.text(textw,texth,textstr,transform=ax.transAxes)

        return h1, ax

    def hist(self,nbins = 20, ylims=None, printstats=True,\
                color='0.5', **kwargs):
        """
        Histogram plot of the errors
        """
        if ylims is None:
            ylims = self.ylims

        bins = np.linspace(ylims[0], ylims[1], nbins)

        textw=0.05
        texth=0.55
        ax = plt.gca()

        h1 = plt.hist(self.error, bins, color=color, **kwargs)
        plt.grid(b=True)
        plt.xlabel('Error [$X_{mod} - X_{obs}$]')
        plt.ylabel('PDF')
        ax.set_xlim(ylims)

        plt.title('StationID: %s'%self.stationid)

        if printstats:
            textstr = 'skill = %6.2f\n$r^2$ = %6.2f\nRMSE = %6.2f\nBias = %6.2f\n'%(\
                self.skill.mean(),self.cc.mean(),self.rmse.mean(),self.bias.mean())
            plt.text(textw,texth,textstr,transform=ax.transAxes)

        return h1, ax

    def qqplot(self, percentiles=[1.,5.,25.,50.,75.,95.,99.],\
                ylims=None, **kwargs):
        """
        Quantile-quantile plot
        """
        idx = operator.and_( ~np.isnan(self.TSmod.y), ~np.isnan(self.TSobs.y) )
        q_mod = np.percentile(self.TSmod.y[idx], percentiles)
        q_obs = np.percentile(self.TSobs.y[idx], percentiles)

        if ylims is None:
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
        #idxobs = ~np.isnan(self.TSobs.y)
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

        #if header:
        #    outstr += "-------------------------------------------------------------- \n"
        #    outstr += "           Mean Model Mean Obs. Std. Mod. Std Obs RMSE   CC    skill \n"
        #    
        #    outstr += "---------- ---------- --------- --------- ------- ------ ----- ------ \n"

        #outstr += " %s [%s]  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  %6.3f  \n"%\
        #    (self.stationid, self.units, self.meanMod,\
        #     self.meanObs, self.stdMod, self.stdObs,\
        #     self.rmse,self.cc,self.skill)



        if header:
            outstr += "|      | Mean Model | Mean Obs. | Std. Mod. | Std. Obs. | RMSE |   CC   | Skill |\n"
            
            outstr += "|------| ---------- | --------- | --------- | ------- | --- | ----- | ------| \n"

        outstr += "| %s [%s] | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (self.stationid,self.units, self.meanMod, self.meanObs, self.stdMod, self.stdObs,\
            self.rmse,self.cc,self.skill)

        if f == None:
            print(outstr)
        else:
            f.write(outstr)

    def printStats2(self,f=None,header=True, stationstr=None):
        """
        Prints the statistics to a markdown language style table
        """
        if 'meanMod' not in self.__dict__:
            self.calcStats()

        outstr=''


        if header:
            outstr += "| Station | Std. Mod. | Std. Obs. | Bias | RMSE |   CC   | Skill |\n"
            
            outstr += "|---------| --------- | --------- | ------- | --- | ----- | ------| \n"

        if stationstr is None:
            stationstr = '%s [%s]'%(self.stationid, self.units)

        outstr += "| %s  | %6.3f | %6.3f |  %6.3f | %6.3f | %6.3f | %6.3f | \n"%\
            (stationstr,  self.stdMod, self.stdObs, self.bias,\
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


class ModVsObsUV(ModVsObs):
    """
    Special validation class for velocity data
    """
    
    def __init__(self, tmod, uv_mod, tobs, uv_obs, **kwargs):
        """
        Pass complex variables (u+iv) instead of real
        """
        
        # Make sure the arrays are complex
        assert uv_mod.dtype == np.dtype('complex128')
        assert uv_obs.dtype == np.dtype('complex128')
        
        ModVsObs.__init__(self, tmod, uv_mod, tobs, uv_obs, **kwargs)
        
        u_hat_mod = uv_mod
        u_hat_obs = uv_obs

        self.speed_mod = np.abs(u_hat_mod)
        self.speed_obs = np.abs(u_hat_obs)

        self.theta_mod = np.angle(u_hat_mod)
        self.theta_obs = np.angle(u_hat_obs)
    
    def plot_speeddirn(self, speedmax=1.5,\
        legend=True, loc='upper right', modcolor='r', obscolor='b', **kwargs):
        
        time = uc.TSmod.t
        
        ax1 = plt.subplot(211)
        plt.plot(time, speed_obs, color='b', linewidth=0.25)
        plt.plot(time, speed_mod, color='r', linewidth=0.25)
        ax1.set_ylim(0,speedmax)
        ax1.set_xticklabels([])
        plt.ylabel('Speed [m/s]')

        if legend:
            plt.legend(('Model','Observed'),loc=loc)


        ax = plt.subplot(212)
        plt.plot(time, theta_obs*180./np.pi, color=obscolor,\
                markersize=0.5, marker='.',**kwargs)
        plt.plot(time, theta_mod*180./np.pi, color=modcolor,\
                markersize=0.5, marker='.',**kwargs)
        
        ax.set_yticks([-180,-90.,0,90, 180])
        ax.set_yticklabels(['W','S','E','N','W'])
        ax.set_ylim(-180,180)
        
        plt.xticks(rotation=17)
        
        return ax1, ax
    
    def plot_pdf_polar(self, speedmax=1.0, cmap='RdYlBu_r'):
        """
        Polar plot of speed-direction joint frequency distribution
        """
        # Create a 2d histogram of speed direction
        abins = np.linspace(-np.pi, np.pi, 90.0)      # 0 to 360 in steps of 360/N.
        sbins = np.linspace(0.0, speedmax, 25) 

        Hmod, xedges, yedges = np.histogram2d(self.theta_mod, self.speed_mod, bins=(abins,sbins), normed=True)
        Hobs, xedges, yedges = np.histogram2d(self.theta_obs, self.speed_obs, bins=(abins,sbins), normed=True)

        #Grid to plot your data on using pcolormesh
        theta = 0.5*abins[1:]+0.5*abins[0:-1]
        r = 0.5*sbins[1:]+0.5*sbins[0:-1]

        # Make sure plot covers the whole circle
        theta[0] = -np.pi
        theta[-1] = np.pi
        #theta, r = np.mgrid[0:2*np.pi:360j, 1:100:50j]

        # Contour levels precentages
        clevs = [0.001, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]
        xlabels = ['E','ENE','N','WNW','W','WSW','S','ESE']

        fig = plt.figure( figsize = (14,6))

        ax1 = plt.subplot(121, projection='polar')
        ax1.set_xticklabels(xlabels)
        C = ax1.contourf(theta, r, Hobs.T, clevs, cmap = cmap)
        #plt.colorbar(C)
        plt.title('StationID: %s - Observed\n'%self.stationid)

        ax2 = plt.subplot(122, projection='polar')
        ax2.set_xticklabels(xlabels)
        C = ax2.contourf(theta, r, Hmod.T, clevs, cmap = cmap)
        #plt.colorbar(C)
        plt.title('StationID: %s - Modelled\n'%self.stationid)

        # Insert a colorbar
        axcb = plt.axes([0.43,0.1, 0.16,0.05])
        plt.colorbar(C, cax=axcb, orientation='horizontal', format='%3.1f')
        axcb.set_title('[% frequency]')

        plt.tight_layout()
        
        return fig, ax1, ax2, axcb


##########
# User functions
##########
def load_netcdf(ncfile, group):
    """
    Load a ModVsObs object previously saved in a netcdf file
    """
    # Load the group into a ModVsObs object
    ds = xray.open_dataset(ncfile, group=group)

    # work out the varname
    varnames = list(ds.data_vars.keys())

    for vv in varnames:
        if 'mod' in vv:
            varname = vv.strip('_mod')

    # Load the two variables (as Pandas objects)
    TSobs = ds['%s_obs'%varname].to_pandas()
    TSmod = ds['%s_mod'%varname].to_pandas()

    # Load the attributes
    attrs = ds['%s_obs'%varname].attrs

    # Get the depth
    try:
        Z = TSobs.columns.values
    except:
        Z = 0.

    # Convert to a ModVsObs object
    # Put the data into a ModVsObs object (model first then observed)
    return ModVsObs(\
    	    TSmod.index.to_pydatetime(),\
	    #TSmod.index.values,\
            TSmod.values,\
            TSobs.index.to_pydatetime(),\
	    #TSobs.index.values,\
            TSobs.values,\
            varname=varname,\
            long_name=attrs['long_name'], \
            units=attrs['units'], \
            stationid=group,\
            Z=Z,\
        )
 
