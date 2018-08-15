
"""
Suite of tools for creating and interacting with netcdf files.

Also includes functions for writing to, and querying an sql database
containing the netcdf metadata.

See the examples at the end of this script for help.

Matt Rayson
Stanford University
August 2012
"""

from netCDF4 import Dataset, num2date
#import shapefile
import numpy as np
#try:
#    from pyspatialite import dbapi2 as db
#    print('Warning: pyspatialite not installed - reverting to sqlite3 library...')
#
#except:
import sqlite3 as db

from datetime import datetime
import matplotlib.pyplot as plt

import xray
import yaml
import os
import soda

import pdb

_path = os.path.abspath(soda.__file__)
_dirpath = os.path.dirname(_path)

# Load the netcdf metadata for all of the routines
with open('%s/dataio/ncmetadata.yaml'%_dirpath, 'r') as f:
    ncmeta = yaml.load(f)

mysqldir = '/home/suntans/.conda/envs/soda/lib/'

###################
# xray routines
###################
def dict_toxray(data, ds={}, **kwargs):
    """
    Converts a dictionary with keys as variable names to an
    xray.Dataset object

    The dictionary keys should correspond with variable names in 
    ncmetadata.yaml

    **kwargs are passed directly to xray.DataArray()
    """
    
    for vv in list(data.keys()):
        
        if vv in ncmeta:
            attrs = ncmeta[vv]['attributes']
        else:
            print('Warning variable: %s not in ncmetadata.yaml. Dataset will have no attrs')
            attrs = {}

        da = xray.DataArray(data[vv], attrs = attrs, **kwargs)

        ds.update({vv:da})

    return xray.Dataset(ds)


def load_sql_ncstation(dbfile, station_name, varname, otherquery=None, query_only=False):

    """
    Load netcdf station data referenced in an sql database file
    
    Inputs:
        dbfile - location of database file
        station_name - StationName in database
        varname - variable name e.g. 'waterlevel', 'discharge', 'salinity'
        otherquery - (optional) string with other query variables
                e.g 'time_start > "2007-01-01 00:00:00" and time_end < "2008-01-01 00:00:00"'
        
           
    Returns:
        an xray.DataArray object
        -1 on error
            
    """
    
    outvar = ['NetCDF_Filename',\
        'NetCDF_GroupID',\
        'StationName',\
        'StationID',\
        'X',\
        'Y',\
        'height_start',\
        'Variable_Name',\
        'time_start',\
        'time_end',\
        ]

    tablename = 'observations'
    if not otherquery is None:
        #condition = 'Variable_Name = "%s" and StationName LIKE "%%%s%%" and %s'\
        #    %(varname, station_name, otherquery)
        condition = 'LOWER(Variable_Name) LIKE LOWER("%s") and StationName LIKE "%%%s%%" and %s'\
            %(varname, station_name, otherquery)
	
        print(condition)

    else:
        #condition = 'Variable_Name = "%s" and StationName LIKE "%%%s%%"'\
        condition = 'LOWER(Variable_Name) LIKE LOWER("%s") and StationName LIKE "%%%s%%"'\
            %(varname, station_name)
    
    # Query the database
    print('Querying database...')
    print(condition)
    query = returnQuery(dbfile,outvar,tablename,condition)

    if query_only:
        return query

    # Loop through and extract the variable datasets from each of the files
    data = []
    ii = 0
    for  ncfile, ncgroup, ncvarname in \
                zip(query['NetCDF_Filename'], query['NetCDF_GroupID'], query['Variable_Name']):
        print(ncfile, ncgroup)

        nc = xray.open_dataset(ncfile, group=ncgroup)

        ncdata = nc[ncvarname]

        # Make sure the metadata is in the attributes
        ncdata.attrs.update({'X':query['X'][ii],\
                'Y':query['Y'][ii],\
                'Z':query['height_start'][ii],\
                'time_start':query['time_start'][ii],\
                'time_end':query['time_end'][ii],\
                'StationName':query['StationName'][ii],\
                'StationID':query['StationID'][ii],\
        })

        data.append(ncdata)

        ii += 1

    return data

###################
# Old routines
###################
def writePointData2Netcdf(ncfile,data,globalatts):
    """ 
    Function for writing point/observation data to a grouped netcdf file.
    Each variable is written to a separate group to allow for variable time and/or
    spatial coordinates.
    
    Inputs: 
        ncfile - name of the output netcdf file (string)
        data - list of dictionaries with netcdf data and attribute info
            (see noaatools.py --> extractIOOS for an example of this array)
        globalatts - dictionary with global attributes

    """    
    
    # Convert the dictionary array to a grouped netcdf file
    print('################################################')
    print(' Writing to file: %s...' % ncfile)
    nc = Dataset(ncfile, 'w', format='NETCDF4')
    # Write the global attributes
    for gg in list(globalatts.keys()):
        nc.setncattr(gg,globalatts[gg])
    
    # Each listing in the dictionary is treated as a separate group
    ctr=-1
    for dd in data:
        # Each variable in the group
        for vv in dd:
            # Create a group 
            ctr = ctr + 1
            groupID = 'groupID_%04d' % ctr
            grp = nc.createGroup(groupID)
            
            # Work out the coordinates and create the dimensions
            for cc in dd[vv]['coords']:
                dimname = cc['Name']
                dimlength = np.size(cc['Value']) 
                grp.createDimension(dimname,dimlength)
                print(dimname, dimlength)
                
                # Create the coordinate variables
                tmpvar=grp.createVariable(cc['Name'],'f8',(dimname,))
                tmpvar[:] = cc['Value']

                
                # Create the attributes
                for aa in list(cc.keys()):
                    if aa !='Name' and aa !='Value':
                        tmpvar.setncattr(aa,cc[aa]) 
            # Now create the varible and attribute data
             
            # The dimension info is stored in the coordinates attribute
            coordList = [str(x) for x in dd[vv]['coordinates'].split(', ')]
            tmpvar = grp.createVariable(vv,'f8',(coordList))
            # Write the data
            print(vv, np.size(dd[vv]['Data']), coordList)
            tmpvar[:] = dd[vv]['Data']
            # Write the attriute data
            for aa in list(dd[vv].keys()):
                if aa !='Data' and aa !='coords':
                    tmpvar.setncattr(aa,dd[vv][aa]) 
    
    nc.close()
    print('Completed writing file.')    
    print('################################################')        
    return

def pointNC2shp(ncfile,shpfile):

    """ Create a shapefile of a netcdf file metadata """
    
    w = shapefile.Writer(shapefile.POINT)
    w.field('long_name')
    w.field('StationName')
    w.field('StationID')
    
    # open the netcdf file
    nc = Dataset(ncfile,'r', format='NETCDF4')
    
    # Loop through the groups
    for grp in nc.groups:
        # Loop through the variables
        for vv in nc.groups[grp].variables:
            # Find the variables that have a coordinates
            allatts = nc.groups[grp].variables[vv].ncattrs()
            if 'coordinates' in allatts:
                # Use this variable
                lon = nc.groups[grp].variables['longitude'][:]
                lat = nc.groups[grp].variables['latitude'][:]
                if 'long_name' in allatts:
                    long_name = nc.groups[grp].variables[vv].long_name
                else:
                    long_name='none'
                if 'StationID' in allatts:
                    StationID = nc.groups[grp].variables[vv].StationID
                else:
                    StationID='none'
                if 'StationName' in allatts:
                    StationName = nc.groups[grp].variables[vv].StationName
                else:
                    StationName='none'
                
                # Write the station data
                w.point(lon,lat)
                w.record(long_name,StationName,StationID)
    
    w.save(shpfile)
    print('NetCDF metadata written to shapefile: %s'%shpfile)
    return

def db2shp(dbfile,shpfile):
    """Converts a database file to a shape file"""
    
    #Initialise the shapefile    
    w = shapefile.Writer(shapefile.POINT)
    w.field('long_name')
    w.field('StationName')
    w.field('StationID')
    w.field('start_time')
    w.field('end_time')
    
    # Connect to the database
    conn = db.connect(dbfile)
    c = conn.cursor()

    c.execute('select "X","Y","Variable_Name", "StationName","StationID","time_start","time_end" from observations')
    #c.execute('select "X" from observations')
    
    for row in c:
        #x.append(float(row[0]))
        #y.append(float(row[1]))
        #long_name.append(row[2])
        #StationName.append(row[3])
        #StationID.append(row[4])
        # Write the station data
        w.point(row[0],row[1])
        w.record(row[2],row[3],row[4],row[5],row[6])

    w.save(shpfile)
        
    c.close()
    return    
    
def createObsDB(dbfile):
    """ Create a database for storing observational netcdf metadata"""
    
    conn = db.connect(dbfile)

    # Add spatialite module
    #   See this banter here:
    #       https://groups.google.com/forum/#!topic/spatialite-users/o0jUwMUqx_g
    conn.enable_load_extension(True) 
    conn.execute("SELECT load_extension('mod_spatialite')") 
    


    c = conn.cursor()
    # Create table
#    tablefields = {'NetCDF_Filename':'text','NetCDF_GroupID':'text','Variable_Name':'text',\
#    'longitude':'real','latitude':'real','lon_start':'real',\
#    'lon_end':'real','lat_start':'real','lat_end':'real',\
#    'time_start':'text','time_end':'text','height_start':'real',\
#    'height_end':'real','StationName':'text','StationID':'text'}
    tablefields = [\
     # 'id',\
     'NetCDF_Filename','NetCDF_GroupID','Variable_Name',\
     'X','Y',\
    'lon_start','lon_end','lat_start','lat_end','time_start','time_end','height_start',\
    'height_end','StationName','StationID',\
        #'GEOMETRY',\
    ]
    fieldtype = [
    #'INTEGER NOT NULL PRIMARY KEY',\
    'text','text','text',\
    'real','real',\
    'real','real','real','real',\
    'text','text','real','real','text','text',\
        #'text',\
    ]

    tablename = 'observations'

    
    ######
    # PySpatialite specific info
    ######
    # initializing Spatial MetaData
    # using v.2.4.0 this will automatically create
    # GEOMETRY_COLUMNS and SPATIAL_REF_SYS
    sql = 'SELECT InitSpatialMetadata()'
    c.execute(sql)

    # creating a POINT table
    #sql = 'CREATE TABLE test_pt ('
    #sql += 'id INTEGER NOT NULL PRIMARY KEY,'
    #sql += 'name TEXT NOT NULL)'
    #c.execute(sql)
    ###########
    # Original info
    ###########
    
    # Create a string to create the table
    tablestr='('
    for ff,tt in zip(tablefields,fieldtype):
        tablestr += ff+' '+tt+','
    tablestr = tablestr[:-1] + ')'
    
    print(tablestr)
    createstr = 'CREATE TABLE %s %s' % (tablename,tablestr)
    c.execute(createstr)

    # creating a POINT Geometry column
    sql = "SELECT AddGeometryColumn('%s',"%tablename
    sql += "'GEOMETRY', 4326, 'POINT', 'XY')"
    c.execute(sql)
    
    c.close()
    return  
    
def netcdfObs2DB(ncfile, dbfile, nctype=1):
    """
    Extract the relevant metadata from a netcdf-4 file with groups

    nctype = 1 : My original format
    nctype = 2 : RPS-WEL converted format
    ...
    """
    def get_meta_type1(nc, grp ,vv):
        # Find the variables that have coordinates
        allatts = nc.groups[grp].variables[vv].ncattrs()

        if 'coordinates' in allatts or 'units' in allatts:
            # Use this variable
            lon = nc.groups[grp].variables['longitude'][:]
            lat = nc.groups[grp].variables['latitude'][:]
            if 'long_name' in allatts:
                long_name = nc.groups[grp].variables[vv].long_name
            else:
                long_name='none'
            if 'StationID' in allatts:
                StationID = nc.groups[grp].variables[vv].StationID
            else:
                StationID='none'
            if 'StationName' in allatts:
                StationName = nc.groups[grp].variables[vv].StationName
            else:
                StationName='none'
                
            try:
                ele = nc.groups[grp].variables['elevation'][:]
            except:
                ele=[0.0]
                
            times = nc.groups[grp].variables['time']
            dates = num2date([times[0],times[-1]],units=times.units)

        return lon, lat, long_name, StationID, StationName, ele, dates

    def get_meta_type2(nc, grp ,vv):
        """
        Processed RPS files
        """

        if vv in ['time','Longitude','Latitude','DepthHeight']:
            return None, None, None, None, None, None, None

        # Use this variable
        try:
            lon = nc.groups[grp].variables['Longitude'][:]
            lat = nc.groups[grp].variables['Latitude'][:]
        except:
            lon = 0.
            lat = 0.
        try:
            long_name = nc.groups[grp].variables[vv].long_name
        except:
            long_name = ''

        StationID='none'
        StationName = nc.groups[grp].stationname
            
        ele = 0.
        if 'DepthHeight' in list(nc.groups[grp].variables.keys()):
             ele += nc.groups[grp].variables['DepthHeight'][:]

        if hasattr(nc.groups[grp].variables[vv], 'height'):
             ele += nc.groups[grp].variables[vv].height

        #try:
        #    ele = nc.groups[grp].variables[vv].height + \
        #        nc.groups[grp].variables['DepthHeight']
        #except:
        #    pdb.set_trace()
        #    ele = 0.0
            
        times = nc.groups[grp].variables['time']
        dates = num2date([times[0],times[-1]],units=times.units)

        return lon, lat, long_name, StationID, StationName, ele, dates

    def get_meta_type3(nc, grp ,vv):
        """
        Processed xray data

        This is probably the ideal way to store metadata for single
        point time series data
        """

        if vv in ['time','Longitude','Latitude','DepthHeight']:
            return None, None, None, None, None, None, None

        # Use this variable
        try:
            lon = nc.groups[grp].getncattr('Longitude')
            lat = nc.groups[grp].getncattr('Latitude')
        except:
            lon = 0.
            lat = 0.

        try:
            long_name = nc.groups[grp].variables[vv].long_name
        except:
            long_name = vv

        StationID=nc.groups[grp].StationID
        StationName = nc.groups[grp].StationName
            
        try:
            ele = - nc.groups[grp].getncattr('Depth') + \
                nc.groups[grp].getncattr('InstrumentDepth')
        except:
            ele = 0.0
            
        try:
            times = nc.groups[grp].variables['time']
            dates = num2date([times[0],times[-1]],units=times.units)
        except:
            times = nc.groups[grp].getncattr('Time')
            dates = [datetime.strptime(times, '%Y-%m-%d %H:%M:%S')]

        return lon, lat, long_name, StationID, StationName, ele, dates
 
    def get_meta_type4(nc, grp ,vv):
        """
        Processed IMOS data

        """

        if vv in ['TIME','LONGITUDE','LATITUDE','NOMINAL_DEPTH']:
            return None, None, None, None, None, None, None

        # Use this variable
        lon = nc.groups[grp].variables['LONGITUDE'][:]
        lat = nc.groups[grp].variables['LATITUDE'][:]
        long_name = nc.groups[grp].variables[vv].long_name

        StationID=nc.groups[grp].site_code 
        StationName = nc.groups[grp].stationname

        times = nc.groups[grp].variables['TIME']
        dates = num2date([times[0],times[-1]],units=times.units)
            
        try:
            ele = nc.groups[grp].variables['NOMINAL_DEPTH'][:]
        except:
            ele = nc.groups[grp].instrument_nominal_depth
            
        return lon, lat, long_name, StationID, StationName, ele, dates
 
 
             
    # Write metadata to the sql database
    
    # Open the database
    conn = db.connect(dbfile)

    # Load the spatialite features
    conn.enable_load_extension(True) 
    conn.execute("SELECT load_extension('%s/mod_spatialite')"%mysqldir) 

    c = conn.cursor()
    tablename = 'observations'
    
    # open the netcdf file
    nc = Dataset(ncfile,'r', format='NETCDF4')
    
    # Loop through the groupsind = returnDictInd(ncdict,vv,sta)
    for grp in nc.groups:
        # Loop through the variables
        for vv in nc.groups[grp].variables:

            if vv in ['time','longitude','latitude','elevation']:
                continue
               
            write = True

            if nctype==1:
                lon, lat, long_name, StationID, StationName, ele, dates = \
                    get_meta_type1(nc, grp ,vv)

   
            elif nctype == 2:
                print(grp, vv)
                lon, lat, long_name, StationID, StationName, ele, dates = \
                    get_meta_type2(nc, grp ,vv)

            elif nctype == 3:
                print(grp, vv)
                lon, lat, long_name, StationID, StationName, ele, dates = \
                    get_meta_type3(nc, grp ,vv)

            elif nctype == 4:
                print(grp, vv)
                lon, lat, long_name, StationID, StationName, ele, dates = \
                    get_meta_type4(nc, grp ,vv)

            if lon is None:
                write = False

            else:
                geom = "GeomFromText('POINT("
                geom += "%f " % (lon)
                geom += "%f" % (lat)
                geom += ")', 4326)"

                dbtuple = (ncfile, grp, vv,\
                    lon, lat,\
                    lon, lon,\
                    lat, lat, dates[0], dates[-1],\
                    ele,ele,StationName,StationID,\
                    geom,\
                    #'POINT',\
                    )
 
    
            if write:
                # Create the tuple to insert into the database
                dbstr = '("%s", "%s", "%s",'
                dbstr += '%4.6f, %4.6f,'
                dbstr += '%4.6f, %4.6f, %4.6f, %4.6f, "%s", "%s",'
                dbstr += '%4.6f, %4.6f, "%s", "%s", %s)'
                dbstr=dbstr%dbtuple

                # Insert into the db
                #execstr = 'INSERT INTO %s VALUES %s'%(tablename,dbstr)
                #print execstr
                c.execute('INSERT INTO %s VALUES %s'%(tablename,dbstr))
                # Save (commit) the changes
                conn.commit()
                
                
    # Insert a row of data
    #c.execute("INSERT INTO stocks VALUES ('2006-01-05','BUY','RHAT',100,35.14)")
    # Close everything                
    nc.close()
    c.close()
    return

def returnGroup(ncfile,grpid):
    """ Return point data from a grouped netcdf file into a dictionary
    For use with query builder
    **Note that the time variable is returned as an array of datetime objects (see datetime module)
    """
   
    nc = Dataset(ncfile,'r', format='NETCDF4')

    output = {}
    for vv in nc.groups[grpid].variables:
        
        if vv == 'time':
            times=nc.groups[grpid].variables[vv]
            dates = num2date(times[:],units=times.units)
            output.update({vv:dates})
        else:
            #output.update({vv: np.ravel(nc.groups[grpid].variables[vv][:])})
            output.update({vv:nc.groups[grpid].variables[vv][:]})
    
    nc.close()
    return output

def returnGroupFast(ncfile,grpid,nc):
    """ Return point data from a grouped netcdf file into a dictionary
    For use with query builder
    **Note that the time variable is returned in its netcdf units (faster)
    """
   
    #nc = Dataset(ncfile,'r', format='NETCDF4')

    output = {}
    for vv in nc.groups[grpid].variables:
        if vv == 'time':
            times=nc.groups[grpid].variables[vv]
            dates = num2date(times[:],units=times.units)
            output.update({vv:dates})
        else:
            #output.update({vv: np.ravel(nc.groups[grpid].variables[vv][:])})
            output.update({vv:nc.groups[grpid].variables[vv][:]})
    
    #nc.close()
    return output
    
def returnQuery(dbfile,outvar,tablename,condition):
    """Returns a dictionary with the fields specified in a query
    
    Example condition:
        'Variable_Name = "RH" and start_time >= "2011-01-01 00:00:00"'
    
    """    
    
    # Open the database
    conn = db.connect(dbfile)
    c = conn.cursor()
    
    querystr = 'SELECT %s FROM %s WHERE %s'%(', '.join(outvar),tablename,condition)
    #print querystr
    query = c.execute(querystr)
    
    output = {}
    for vv in outvar:
        output.update({vv:[]})
    
    
    for row in query:
        k=-1
        for vv in outvar:
            k+=1
            output[vv].append(row[k])
    c.close()
    return output

def queryNC(dbfile,outvar,tablename,condition,fastmode=False):
    """ Main function for extracting queried data for netcdf files
    
    Example inputs:
    dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
    outvar = ['NetCDF_Filename','NetCDF_GroupID']
    tablename = 'observations'
    condition = 'Variable_Name = "RH"'
    """    
    # Get the query
    query = returnQuery(dbfile,outvar,tablename,condition)
    
    data=[]
    ii=-1
    for ff,gg in zip(query['NetCDF_Filename'],query['NetCDF_GroupID']):
        ii+=1
        #print 'Extracting data from: %s, Group: %s...'%(ff,gg) 
        if fastmode:
            if ii == 0:
                nc = Dataset(ff,'r', format='NETCDF4')
                ffold = ff
            if ffold != ff:
                nc.close()
                nc = Dataset(ff,'r', format='NETCDF4')
                ffold = ff
            data.append(returnGroupFast(ff,gg,nc))
        else:
            data.append(returnGroup(ff,gg))
    if fastmode:
        nc.close()
    return data, query

def createObsDict(varname,longname,units,data,time,latitude,longitude,height,stationid,stationname,ncdict=[] ):
    """
    Create the list of dictionaries expected by writePointData2Netcdf
    
    All inputs are lists of arrays except: varname, longname, units should be strings
    
    Can append to a previously created dictionary by setting ncdict
    """
    
    # Initialise the output dictionary
    for ID,nn,lat,lon,hh,tt,dd in zip(stationid,stationname,latitude,longitude,height,time,data):

        coords=[{'Name':'longitude','Value':lon,'units':'degrees East'},\
            {'Name':'latitude','Value':lat,'units':'degrees North'},\
            {'Name':'elevation','Value':hh,'units':'metres','positive':'up'},\
            {'Name':'time','Value':tt,'units':'minutes since 1970-01-01 00:00:00'}]
        # Can't put elevation in the "coordinates" list as it expands the array
        attribs = {'StationID':ID,'StationName':nn,'Data':dd,\
            'coordinates':'time, elevation, latitude, longitude','long_name':longname,\
            'units':units,'coords':coords} 
        ncdict.append({varname:attribs})
            
    return ncdict

# Examples

####
# Example 1) Writing metadata to a database
####
def runExample1():
    #ncfile = 'C:/Projects/GOMGalveston/DATA/Winds/NCDCNWS_AirObs_2011.nc'
    ncfile = 'C:/Projects/GOMGalveston/DATA/Ocean/USIOOS_OceanObs_20102011.nc'
    dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
    #createObsDB(dbfile)
    netcdfObs2DB(ncfile,dbfile)
    
    
####
# Example 2) Return data from all stations in a database with Variable_Name = varname
####
def runExample2():

    import matplotlib.dates as dates

    dbfile = 'C:/Projects/GOMGalveston/DATA/GalvestonObs.db'
    outvar = ['NetCDF_Filename','NetCDF_GroupID','StationName']
    tablename = 'observations'
    varname = 'waterlevel'
    condition = 'Variable_Name = "%s"' % varname
    
    data, query = queryNC(dbfile,outvar,tablename,condition)
    
    # Plot the results in one figure
    datemin = datetime(2011,6,1)
    datemax = datetime(2011,7,1)
    ylim = [27,35] # temp
    ylim = [0, 100] 
    fig = plt.figure(figsize=(8,12))
    plt.hold(True)
    k=0
    for dd in data:
        # convert the time for plotting
        t = dates.date2num(dd['time'])
        k+=1
        ax = plt.subplot(len(data),1,k)
        plt.plot(t,dd[varname],'b')
        plt.title('%s at %s'%(varname,query['StationName'][k-1]))
        # Format the x-ticks
        ax.set_xlim(dates.date2num(datemin), dates.date2num(datemax))
        ax.set_ylim(ylim[0],ylim[1])
        ax.grid(True)
        ax.xaxis.set_major_formatter(dates.DateFormatter('%d%b%Y'))
      
    fig.autofmt_xdate() 
    plt.show()
              
