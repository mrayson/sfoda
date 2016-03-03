"""
Tools for downloading specfific ocean/atmosphere/climate
datasets from an opendap server
"""

from mythredds import GetDAP, Dataset, MFncdap
from datetime import datetime,timedelta

import pandas as pd
import numpy as np

###
# Dictionary containing opendap server specific
# information for each model
###
metoceandict = {\
    'HYCOM':{\
        'ncurl':'http://tds.hycom.org/thredds/dodsC/glb_analysis',\
        'type':'ocean',\
        'u':'u',\
        'v':'v',\
        'temp':'temperature',\
        'salt':'salinity',\
        'ssh':'ssh',\
    },\
    'HYCOM_REANALYSIS':{\
        'ncurl':'http://tds.hycom.org/thredds/dodsC/GLBu0.08/expt_19.1/3hrly',\
        'type':'ocean',\
        'u':'water_u',\
        'v':'water_v',\
        'temp':'water_temp',\
        'salt':'salinity',\
        'ssh':'surf_el',\
    },\
    'BRAN_3p5':{\
        'ncurl':[],\
        'multifile':True,\
        'type':'ocean',\
        'u':'u',\
        'v':'v',\
        'temp':'temp',\
        'salt':'salt',\
        'ssh':'eta_t',\
        'u_file':'u',\
        'v_file':'v',\
        'temp_file':'temp',\
        'salt_file':'salt',\
        'ssh_file':'eta',\
    },\
    'GFS':{\
        'ncurl':[],\
        'type':'atmosphere',\
        'multifile':True,\
        'uwind':'ugrd10m',\
        'vwind':'vgrd10m',\
        'tair':'tmp2m',\
        'pair':'prmslmsl',\
        'rh':'rh2m',\
        'cloud':'tcdcclm',\
        'rain':'pratesfc',\
    },\
    'CFSR_1HR':{\
        'ncurl':[],\
        'type':'atmosphere',\
        'multifile':True,\
        'uwind':'U-component_of_wind',\
        'vwind':'V-component_of_wind',\
        'tair':'Temperature',\
        'pair':'Pressure_reduced_to_MSL',\
        'rh':None,\
        'cloud':None,\
        'rain':'Precipitation_rate',\
        'dlwr':'Downward_Long-Wave_Rad_Flux',\
        'dswr':'Downward_Short-Wave_Rad_Flux',\
        'sh':'Specific_humidity',\
        'uwind_file':'wnd10m',\
        'vwind_file':'wnd10m',\
        'tair_file':'tmp2m',\
        'pair_file':'prmsl',\
        'rh':None,\
        'cloud':None,\
        'rain_file':'prate',\
        'dlwr_file':'dlwsfc',\
        'dswr_file':'dswsfc',\
        'sh_file':'q2m',\
    },\

} # End of dictionary

#############################
# Dataset specific classes 
#############################

# These are for doing thigs like retrieving file names,
# setting coordinate projections etc
class GFSFiles:
    """
    Class that returns the filenames of global forecast system output
    for a given time range
    """

    resolution = 4 # 3 = one degree, 4 = 0.5 degree
    #dt = 3 # time interval between files in hours
    #baseurl =\
        #'http://nomads.ncdc.noaa.gov/thredds/dodsC/gfs-%03d/%s/%s/gfs_%d_%s_0000_%03d.grb2'
    dt = 24 
    baseurl =\
        'http://nomads.ncep.noaa.gov:9090/dods/gfs_0p50/gfs%s/gfs_0p50_00z'

    def __init__self(**kwargs):
        self.__dict__.update(**kwargs)

    def __call__(self,trange):

        self.basetime = datetime(1900,1,1)

        trange = [datetime.strptime(trange[0],'%Y%m%d.%H%M%S'),\
            datetime.strptime(trange[1],'%Y%m%d.%H%M%S')]

        # Generate a list of time variables
        dt = timedelta(hours=self.dt)
        t1 = trange[0]
        time = [t1]
        while t1 <= trange[1]:
            t1+=dt
            time.append(t1)
            print t1

        return time,[self.generate_url(tt) for tt in time]

    def generate_url(self,time):
        """
        Generate a url for a given time
        """
        def _gen_url(yymmdd,yyyymm,hours):
            #return self.baseurl%(self.resolution,\
            #    yyyymm,yymmdd,self.resolution,\
            #    yymmdd,hours)
            return self.baseurl%(yymmdd)


        yymmdd = datetime.strftime(time,'%Y%m%d')
        basetime = datetime.strptime(yymmdd,'%Y%m%d')

        # Generate the string
        yyyymm = datetime.strftime(time,'%Y%m')
        hours = (time-basetime).total_seconds()/3600

        url = _gen_url(yymmdd,yyyymm,hours)

        # Check if the url exists
        if not basetime == self.basetime:
            print 'Checking if url exists...\n\t%s'%url
            try:
                # Update to a new data
                #f = urllib2.urlopen('%s.html'%url)
                nc = Dataset(url)
                self.basetime = basetime
                print 'yes'
                nc.close()
                return url
            except:
                print 'File does not exist - we are in the forecast\
                    stage...(%s)'%(yymmdd)
                # Generate a string from the old basetime 
                yymmdd = datetime.strftime(self.basetime,'%Y%m%d')
                yyyymm = datetime.strftime(self.basetime,'%Y%m')
                hours = (time-self.basetime).total_seconds()/3600
                url = _gen_url(yymmdd,yyyymm,hours)
                return url

class CFSR_1hr(object):
    """
    Class for returning the list of CFSR files for a given time 
    range and variable list
    """

    baseurl = 'http://nomads.ncdc.noaa.gov/thredds/dodsC/cfsr1hr/'
    dt = 1.
    dt_units = 'H'

    time_min = '1979-01-01'
    time_max = '2009-12-31'

    def __init__(self, trange, tdsdict, **kwargs):
        self.__dict__.update(**kwargs)

        self.tdsdict = tdsdict

        # generate a list of time variables
        self.t0 = pd.datetime.strptime(trange[0], '%Y%m%d.%H%M%S')
        self.t1 = pd.datetime.strptime(trange[1], '%Y%m%d.%H%M%S')

        time = pd.date_range(self.t0, self.t1,\
               freq = '%d%s'%(self.dt,self.dt_units))

        self.time = np.array([tout.to_datetime() for tout in time])

        #self.timestr = [pd.datetime.strptime(tt, '%Y%m%d.%H%M%S') for tt in self.time]

        # Create a dictionary with the file names for the first variables only
        fstr = self.tdsdict['%s_file'%'uwind']
        filelist = [self.generate_url(self.baseurl, tt, fstr) for tt in self.time]

        self.ncfilelist = np.unique(filelist)
        self.tdsdict['ncurl'] = self.ncfilelist

        # Create an MFncdap object for time/variable lookup
        self.MF = MFncdap(self.ncfilelist)

    def __call__(self, time, var=None):
        """
        Return the same output as a call to MFncdap:
            filename, tstep
        """
        tind, fnames,tslice = self.MF(time)

        # Replace the file name with the correct one for my variable
        fstr = self.tdsdict['%s_file'%'uwind']
        fstrnew = self.tdsdict['%s_file'%var]
        fnamenew = [fname.replace(fstr, fstrnew) for fname in fnames]

        # Replace the keys in the dictionary
        for ff in tslice.keys():
            new_key = ff.replace(fstr, fstrnew)
            tslice[new_key] = tslice.pop(ff)


        return tind, fnamenew, tslice
        #return filename, tslice, vname

    def generate_url(self,baseurl, time, varname):
        """
        Create the url string
        """
        timestr = pd.datetime.strftime(time,'%Y%m')
        return '%s%s/%s.gdas.%s.grb2'%(baseurl,timestr,varname,timestr)

    def get_filename_only(self, var=None):
        """
        Returns the first file only
        """
        fname =  self.ncfilelist[0]
        if not var == None:
            fstr = self.tdsdict['%s_file'%'uwind']
            fstrnew = self.tdsdict['%s_file'%var]
            fnamenew = fname.replace(fstr, fstrnew)
            
        else:
            fnamenew = fname

        return fnamenew

class GetDAP_BRAN(GetDAP):
    """
    Special class for dealing with Bran data
    """
    
    def __init__(self,**kwargs):
        GetDAP.__init__(self, **kwargs)

    def get_coord_names(self,varname):
        """
        Just use the dimension names
        """
        timecoord,zcoord, ycoord,xcoord =  self._nc.variables[varname].dimensions
        return timecoord, xcoord, ycoord,zcoord


##############################
# Main functions for calling
##############################
def get_metocean_dap(xrange,yrange,zrange,trange,outfile,\
        gridfile=None,name='HYCOM',**kwargs):
    """
    Extracts any met-ocean model that doesn't have server specific needs

    Inputs:
        name is dataset name in metoceandict (HYCOM by default)
        xrange/yrange/zrange: lists with i.e. [min(x),max(x)]

        trange: lists with start and end time in string format i.e.
            ['20101201.000000','20101202.000000']
    """
    oceandict=metoceandict[name]

    for key in kwargs.keys():
        oceandict[key]=kwargs[key]

    # Construct the dap class
    TDS = GetDAP(gridfile=gridfile,**oceandict)
    # Call the object
    TDS(xrange,yrange,trange,zrange=zrange,outfile=outfile)

    print 'Done.'
    return TDS


def get_metocean_local(ncfile,varname,name='HYCOM',TDS=None,\
        xrange=None,yrange=None,zrange=None,trange=None):
    """
    Retrieves variable data from a local file

    Use common name for variable i.e. u,v,temp,tair,uwind, etc

    Returns: data, ncobject (containing coordinates, etc)

    This is sort of a long way around but is general
    """
    oceandict=metoceandict[name]

    # 
    ncvar = oceandict[varname]
    oceandict['ncurl']=ncfile
    if oceandict.has_key('multifile'):
        oceandict['multifile']=False

    # Construct the dap class
    if TDS==None:
        if name == 'BRAN_3p5': # This is a total hack but oh well...
            TDS = GetDAP_BRAN(variables=[varname],**oceandict)
        else:
            TDS = GetDAP(variables=[varname],**oceandict)

    ncobj = getattr(TDS,'_%s'%varname)

    # Call the object
    data = ncobj.get_data(ncvar,xrange,yrange,zrange,trange)

    return data, ncobj


def get_gfs_tds(xrange,yrange,zrange,trange,outfile):
    """
    Extract Global forecasting system data
    """
    gfsdict = metoceandict['GFS']

    # Get the file names for the given time range from the class
    gfs = GFSFiles()
    time,files = gfs(trange)

    # Update the dictionary
    gfsdict['ncurl']=files

    # Create the thredds object
    TDS = GetDAP(**gfsdict)
    
    # Call the object
    TDS(xrange,yrange,trange,zrange=zrange,outfile=outfile)

def get_cfsr_tds(xrange, yrange, trange, outfile, outfile_pair):
    """
    Extract Climate forecasting system (CFSR) data
    """
    vars = [
        'uwind',\
        'vwind',\
        'tair',\
        #'pair',\
        'rain',\
        'dlwr',\
        'dswr',\
        'sh',\
        ]

    mydict = metoceandict['CFSR_1HR']
    cfsr = CFSR_1hr(trange, mydict)

    # Create the thredds object
    TDS = GetDAP(vars = vars, MF = cfsr, **mydict)
    # Call the object
    TDS(xrange,yrange,trange,outfile=outfile)

    # Note that pressure is on a separate grid so we will store it separately
    TDS = GetDAP(vars = ['pair'], MF = cfsr, **mydict)
    # Call the object
    TDS(xrange,yrange,trange,outfile=outfile_pair)



