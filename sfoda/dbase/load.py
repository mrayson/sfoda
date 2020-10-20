
import numpy as np
from .netcdfio import queryNC

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


