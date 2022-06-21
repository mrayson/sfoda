# coding: utf-8

"""
Process RDI adcp data using dolfyn (lkilcher.github.io/dolfyn/)
"""

#try:
#    from mydolfyn.adp import api 
#    #from dolfyn.adp import rotate
#    from .transform import Transform, rdi_xyz_enu
#
#except:
#    print('Warning - dolfyn not installed. ADCP io not available.')

#from dolfyn.adp import api 
from dolfyn.io import rdi
#from dolfyn.adp import rotate
from .transform import Transform, rdi_xyz_enu


import xarray as xray
import os
from datetime import datetime, timedelta
from collections import OrderedDict
import numpy as np
from scipy.interpolate import interp1d

import pdb

deg2rad = np.pi / 180.
rad2deg = 180./np.pi

class ADCP_io(object):
    """
    Wrapper object for parsing an RDI file with the "dolfyn" package from github

    Attributes:
    --------
        ds   - xray dataset with the pertinent data
        _data - adcp data object from dolfyn

    Methods:
    --------
        to_netcdf : save xray.Dataset to netcdf

    """

    
    def __init__(self, adcpfile, group=None, rotate=False, \
        mapbins=False, config_update = {}):
        """
        Inputs:
        -----
                adcpfile - RDI binary or a netcdf file
                rotate - [False] option to rotate
                config_update - dictionary with keys to replace in the ADCP config file

        """

        self.adcpfile = adcpfile

        self.rotate = rotate
        self.mapbins = mapbins

        try:
            # Load data directly from a netcdf file
            self.ds = xray.open_dataset(self.adcpfile, group=group, decode_coords=False)

        except:
            # Load the raw binary data
            #self._data = rdi.read_rdi(self.adcpfile)
            #self._data = rdi.adcp_loader(self.adcpfile)
            userdata = rdi.read_userdata(self.adcpfile, None)
            with rdi.adcp_loader(self.adcpfile) as ldr:
                 self._data = ldr.load_data()
            self._data['props'].update(userdata)

            # Override settings in the configuration dictionary
            for kk in list(config_update.keys()):
                self._data.config[kk] = config_update[kk]

            self.positive = self._data.config['orientation']
                
            # Rotate raw data
            if rotate:
                print('Rotating...')
                #self._rotate_velocity()
                #self._rotate_velocity_uhdas(self._data._u, \
                self.u, self.v, self.w, self.errvel,\
                self.u_inst, self.v_inst, self.w_inst = \
                    self._rotate_velocity(self._data.vel.astype('float64').swapaxes(0,1), \
                        self._data.range,\
                        self._data.orient.heading.astype('float64'),\
                        self._data.orient.pitch.astype('float64'),\
                        self._data.orient.roll.astype('float64'),\
                        float(self._data.config['beam_angle']),\
                        self._data.config['beam_pattern'],\
                        self._data.config['orientation'],\
                        self._data.config['coord_sys'])

            self._convert_time()

            self.ds, self.encoding = self._to_xray(self._data)

        ## Rotate netcdf data
        #if rotate:
        #    print 'Rotating...'
        #    #self._rotate_velocity_uhdas(self.ds['beamvel'].values, \
        #    self.ds['u'].values, self.ds['v'].values,\
        #        self.ds['w'].values, self.ds['errvel'].values = \
        #        self._rotate_velocity(self.ds['beamvel'].values, \
        #            self.ds['distance'],\
        #            self.ds['heading'].values,\
        #            self.ds['pitch'].values,\
        #            self.ds['roll'].values,\
        #            self.ds.attrs['beam_angle'],\
        #            self.ds.attrs['beam_pattern'],\
        #            self.ds.attrs['orientation'],\
        #            self.ds.attrs['coord_sys'],)

    def to_netcdf(self, ncfile, mode='w', group=None):
        """
        Save to a netcdf file
        """
        print('Saving data to file %s...'%ncfile)

        if 'encoding' in list(self.__dict__.keys()):
            encoding=self.encoding
        else:
            encoding = None

        self.ds.to_netcdf(ncfile,mode=mode, group=group, encoding=encoding)
        print('Done.')

    def resample_uvw(self, dtavg, raw=False, rotate=False, othervars=[]):
        """
        Resample the velocity data by averaging and subsampling at the average interval


        Inputs:
             dtavg - average interval [seconds]
        Options:
             raw - [False] subsample raw fields
             rotate - [False] if raw, rotate as well. 
                Note that: heading/pitch/roll get averaged BEFORE transformation
             othervars - list of other variables to subsample (not u/v/w)
        """

        def resample(dset, dtavg, axis=0):
            dtstr = '%dS'%(int(dtavg))
            phi = (dset.to_pandas()).resample(dtstr, axis=axis).mean()

            # Use interpolate to fill in nan's
            try:
                return xray.DataArray(phi.interpolate(axis=axis) )
            except:
                return xray.DataArray(phi)

        # Create some global attributes
        attrs = {
                'Name':'Sub-sampled ADCP data',
                'Original File':self.adcpfile,
                'Sub-sample interval':dtavg,
        }
        
        ds = self.ds

        dsnew = xray.Dataset(attrs=attrs)

        if raw:
            dsnew.update({'beamvel' : resample( ds['beamvel'], dtavg, axis=2) })

            #####
            # Rotating after subsampling is a bad idea...
            # Leave it in here for now
            #####
            if rotate:
                # Filter the heading/pitch/roll first
                dsnew.update({'heading' : resample( ds['heading'], dtavg, axis=0) })
                dsnew.update({'pitch' : resample( ds['pitch'], dtavg, axis=0) })
                dsnew.update({'roll' : resample( ds['roll'], dtavg, axis=0) })

                u, v, w, evel, uinst, vinst, vinst = \
                    self._rotate_velocity(dsnew['beamvel'], \
                        dsnew['heading'],\
                        dsnew['pitch'],\
                        dsnew['roll'],\
                        ds.attrs['beam_angle'],\
                        ds.attrs['beam_pattern'],\
                        ds.attrs['orientation'],\
                        ds.attrs['coord_sys'],)

                dsnew.update({'u':u, 'v':v, 'w':w,'errvel':evel})
                
        else:
            dsnew.update({'v' : resample( ds['v'], dtavg, axis=1) })
            dsnew.update({'u' : resample( ds['u'], dtavg, axis=1) })
            dsnew.update({'w' : resample( ds['w'], dtavg, axis=1) })

        for vv in othervars:
             dsnew.update({vv : resample( ds[vv], dtavg, axis=0) })

        return dsnew

    #######
    # QA/QC methods
    #######
    def qaqc_pgood(self, thresh):
        """
        Mask u/v/w based on percentage good information

        Uses minimum of all four beams
        """

        ds = self.ds

        pgood = ds['percent_good'].min(axis=1)

        mask = pgood.values < thresh

        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan

        self.ds.attrs.update({'QAQC Percentage Good Threshold':thresh})

    def qaqc_corr(self, thresh):
        """
        Mask u/v/w based on correlation

        Ensure that the minimum of three beams is above threshold

        Use beams 1-3 for now as these are used for the transformation

        TODO: Use the correlation to indicate which beams to use for transformation
        """
        ds = self.ds
        corr = ds['corr'][:,0:3,:].min(axis=1)

        mask = corr.values < thresh

        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan

        self.ds.attrs.update({'QAQC Correlation Minimum Threshold':thresh})

    def qaqc_echo(self, thresh):
        """
        Mask u/v/w based on echo intensity

        Ensure that all four beams are above threshold
        """
        ds = self.ds
        echo = ds['echo'][:,:,:].min(axis=1)

        mask = echo.values < thresh

        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan

        self.ds.attrs.update({'QAQC Echo Intensity Minimum Threshold':thresh})

    def qaqc_errvel(self, thresh):
        """
        Mask u/v/w based on error velocity

        """

        ds = self.ds

        evel = ds['errvel']

        mask = evel.values > thresh

        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan

        self.ds.attrs.update({'QAQC Error Velocity Maximum Threshold':thresh})


    def qaqc_tilt(self, cutoff_angle):
        """
        Mask u/v/w arrays when the instrument tilt goes too far
        """
        ds = self.ds

        tilt = self._calc_tilt(ds['pitch'].values,\
               ds['roll'].values)

        flag_angle = tilt > cutoff_angle

        # Create a 2d mask array
        nz = ds.distance.shape
        mask = flag_angle[np.newaxis,:].repeat(nz,0)

        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan

        self.ds.attrs.update({'QAQC Maximum Tilt Angle':cutoff_angle})

    def qaqc_depth(self, max_depth, orientation=None, P=None):
        """
        Mask out regions outside of the maximum depth
        """

        ds = self.ds
        if orientation is None:
            orientation = ds.orientation
        else:
            self.ds.attrs['orientation'] = orientation # Update this

        if 'zhat' not in list(ds.keys()):
            self._calc_depth(orientation=orientation, P=P)

        zhat = self.ds.zhat.values
        z = ds.distance.values
        dz = np.abs(z[1] - z[0])


        if orientation == 'down':
            #flag_depth = zhat > max_depth - 0.15*max_depth # too conservative
            mask = zhat > max_depth - 1.5*dz
        else:
            #flag_depth = zhat < 0. + 0.15*max_depth
            mask = zhat < 0. + 1.5*dz


        ds['u'].values[mask] = np.nan
        ds['v'].values[mask] = np.nan
        ds['w'].values[mask] = np.nan


    #######
    # Private methods
    #######
    def _calc_depth(self, orientation=None, P=None):
        """
        Calculate the bin depths as a 2D array using the instrument pressure
        """
        ds = self.ds

        tilt = self._calc_tilt(ds['pitch'].values,\
               ds['roll'].values)

        if orientation is None:
            orientation = ds.orientation

        # Calculate the depth using the TILT and pressure
        z = ds.distance.values
        dz = np.abs(z[1] - z[0])
        zhat = z[:,np.newaxis] * np.cos(tilt*deg2rad)

        # Depth of the instrument
        if P is None:
            P = ds['pressure'].values


        if ds.orientation == 'down':
            zhat = P + zhat
            #flag_depth = zhat > max_depth - 0.15*max_depth # too conservative
            #flag_depth = zhat > max_depth - 1.5*dz
        else:
            zhat = P - zhat
            #flag_depth = zhat < 0. + 0.15*max_depth
            #flag_depth = zhat < 0. + 1.5*dz


        # Update the internal object
        zda = xray.DataArray(zhat,
            dims = ('distance','time',),\
            coords = {'time':ds.time.values,\
                    'distance':ds.distance.values},
            attrs = {\
                    'long_name':'Depth below free-surface',\
                    'units':'m',\
                    'positive':'down',\
                    },
        )

        ds.update({'zhat':zda})


    def _calc_tilt(self, pitch, roll):
        """
        Calculate instrument tilt
        """
        # See the IMOS wiki for recomendations on this:
        # https://github.com/aodn/imos-toolbox/wiki/QCProcedures#adcp-tilt-test---imostiltvelocitysetqc---compulsory
        #TILT = acos(sqrt(1 - sin(ROLL)^2 - sin(PITCH)^2))

        TILT = np.arccos(np.sqrt(1 - np.sin(roll*deg2rad)**2\
                - np.sin(pitch*deg2rad)**2)) * rad2deg

        return TILT

    def _rotate_velocity(self, beamvel, distance,\
        heading_deg, pitch_deg, roll_deg,\
        beam_angle, beam_pattern, orientation, coord_sys):
        """
        Rotate from beam to compass coordinates

        Using the rotation routines from the dolyfn libary

        """
        if not coord_sys == 'beam':
            print('Data collected in %s coordinates - not rotating.'%(coord_sys))
            # error velocity is stored first
            return beamvel[:,0,:], beamvel[:,1,:], beamvel[:,2,:]
            #return beamvel[:,1,:], -beamvel[:,0,:], beamvel[:,2,:]

        # local copy of data for convenience
        #data = self._data

        # Map onto the correct depth-bins
        if self.mapbins:
            beamvelnew = binmap(beamvel, distance, pitch_deg, roll_deg)
        else:
            beamvelnew = beamvel

        # Calculate the rotation matrix
        #isconvex = data.config['beam_pattern'] == 'convex'
        isconvex = beam_pattern == 'convex'
        rotmat = calc_beam_rotmatrix(theta = beam_angle,convex=isconvex)
        print(rotmat)

        #data.config.update({'rotmat':rotmat})

        ## Update the beam
        #data.add_data('beam1vel',data._u[:,0,:])
        #data.add_data('beam2vel',data._u[:,1,:])
        #data.add_data('beam3vel',data._u[:,2,:])
        #data.add_data('beam4vel',data._u[:,3,:])


        # Rotate the beam to instrument coordinates
        #  this populates the variables: u_inst, v_inst
        u_inst, v_inst, w_inst, errvel = beam2inst(beamvelnew, rotmat)

        # Rotate the instrument to earth coordinates
        #self.u, self.v, self.w = rotate.inst2earth(data)

        u, v, w = inst2earth(u_inst, v_inst, w_inst,\
            heading_deg, pitch_deg, roll_deg,\
            orientation=orientation,\
            heading_offset=None, declination=None,\
            fixed_orientation=False)

        #return u_inst, v_inst, w_inst
        return u, v, w, errvel, u_inst, v_inst, w_inst

    def _rotate_velocity_uhdas(self, beamvel,\
        heading_deg, pitch_deg, roll_deg,\
        beam_angle, beam_pattern, orientation):
        """
        Rotate from beam to compass coordinates

        Using the rotation routines from the uhdas/pycurrents library

        """

        # Create the transform objecto
        tr = Transform(angle=beam_angle, geometry=beam_pattern)

        # transform the beam to instrument coordinate
        bvel = beamvel.swapaxes(0,1).swapaxes(0,2)
        xyze = tr.xyz_to_beam(bvel)

        # rotate to east, north, up coordinates
        uvwe = rdi_xyz_enu(xyze,\
                heading_deg, pitch_deg, roll_deg,\
                orientation=orientation)

        return uvwe[:,:,0].T, uvwe[:,:,1].T, uvwe[:,:,2].T

    def _convert_time(self):
        """ Convert the time"""

        # It doesn't appear to be matlab datenum format
        def adcptime2python(matlab_datenum):
            return datetime.fromordinal(int(matlab_datenum)) +\
                timedelta(days=matlab_datenum%1)# - timedelta(days = 366)

        self.time = [adcptime2python(tt) for tt in self._data.mpltime]

    def _to_xray(self, data):
        """
        Create a dictionary of the output variables 
        and convert to an xray Dataset
        """
        #data = self._data

        # dimensions
        adcp_dims = {
            'distance':data.range.squeeze(),
            'time':self.time,
            'beam':list(range(1,5)),
        }

        # Variables
        adcp_vars = {
            'beamvel':
                {'data':data.vel.swapaxes(0,1),
                 'attrs':{'long_name':'Beam velocity',
                     'units':'m/s'},
                 'dims':('distance','beam','time')
                },
            'percent_good':
                {'data':data.signal.prcnt_gd.swapaxes(0,1),
                 'attrs':{'long_name':'Percentage good',
                     'units':''},
                 'dims':('distance','beam','time')
                },
            'echo':
                {'data':data.signal.echo.swapaxes(0,1),
                 'attrs':{'long_name':'Echo intensity',
                     'units':''},
                 'dims':('distance','beam','time')
                },
            'corr':
                {'data':data.signal.corr.swapaxes(0,1),
                 'attrs':{'long_name':'Correlation',
                     'units':''},
                 'dims':('distance','beam','time')
                },
            'u':
                {'data':self.u,
                 'attrs':{'long_name':'Eastward water velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
            'v':
                {'data':self.v,
                 'attrs':{'long_name':'Northward water velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
            'w':
                {'data':self.w,
                 'attrs':{'long_name':'Vertical water velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
            'uinst':
                {'data':self.u_inst,
                 'attrs':{'long_name':'Instrument x velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
            'vinst':
                {'data':self.v_inst,
                 'attrs':{'long_name':'Instrument y velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
            'winst':
                {'data':self.w_inst,
                 'attrs':{'long_name':'Instrument upward velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },
 
            'errvel':
                {'data':self.errvel,
                 'attrs':{'long_name':'Error velocity',
		     'coordinates':'time distance',
                     'units':'m/s'},
                 'dims':('distance','time')
                },

            'pressure':
                {'data':data.depth_m,
                 'attrs':{'long_name':'Pressure',
                     'units':'decibars'},
                 'dims':('time',)
                },
            'temperature':
                {'data':data.env.temperature_C,
                 'attrs':{'long_name':'Water temperature',
                     'units':'degrees C'},
                 'dims':('time',)
                },
            'heading':
                {'data':data.orient.heading,
                 'attrs':{'long_name':'Instrument heading',
                     'units':'degrees'},
                 'dims':('time',)
                },
            'pitch':
                {'data':data.orient.pitch,
                 'attrs':{'long_name':'Instrument pitch',
                     'units':'degrees'},
                 'dims':('time',)
                },
            'roll':
                {'data':data.orient.roll,
                 'attrs':{'long_name':'Instrument roll',
                     'units':'degrees'},
                 'dims':('time',)
                },
        }

        # Create an output dataset
        # Global attrs
        attrs = data.config

        ds = xray.Dataset(attrs=attrs)

        encoding = {}
        for var in list(adcp_vars.keys()):
        #for var in ['beamvel']:
            print('Converting variable: %s...'%var)
            
            coords = OrderedDict()
            for dd in adcp_vars[var]['dims']:
                coords.update({dd:adcp_dims[dd]})
                
            V = xray.DataArray( adcp_vars[var]['data'],\
                dims=adcp_vars[var]['dims'],\
                name=var,\
                attrs = adcp_vars[var]['attrs'],\
                coords = coords
            )

            ds.update({var:V})
            
            encoding.update({var:{'zlib':True,'_FillValue':-999999.}})

        return ds, encoding


#########
# Instrument rotation routines
#
# Adapted from dolfyn library's .../adp/rotate.py
#########

def calc_rotation_terms(theta, convex, degrees):
    if degrees:
        theta = theta * deg2rad
    if convex:
        c = 1.
    else:
        c = -1.

    a = 1 / (2. * np.sin(theta))
    b = 1 / (4. * np.cos(theta))
    d = a / (2. ** 0.5)

    return a, b, c, d

def calc_beam_rotmatrix(theta=20.0, convex=True, degrees=True):
    """Calculate the rotation matrix from beam coordinates to
    instrument head coordinates.

    Parameters
    ----------
    theta : is the angle of the heads (usually 20 or 30 degrees)

    convex : is a flag for convex or concave head configuration.

    degrees : is a flag which specifies whether theta is in degrees
        or radians (default: degrees=True)
    """
    a,b,c,d = calc_rotation_terms(theta, convex, degrees)

    return np.array([[c * a, -c * a, 0, 0],
                     [0, 0, -c * a, c * a],
                     [b, b, b, b],
                     [d, d, -d, -d]])

    # Consistent with bm2dir.m for an upward looking instrument
    #return np.array([[-c * a, c * a, 0, 0],
    #                 [0, 0, -c * a, c * a],
    #                 [-b, -b, -b, -b],
    #                 [-d, -d, d, d]])



def beam2inst(beamvel, rotmat):
    """Rotate velocities from beam to instrument coordinates.
    """
    #if hasattr(adcpo.config, 'rotmat'):
    #    rotmat = adcpo.config.rotmat
    #else:
    #    rotmat = calc_beam_rotmatrix(adcpo.config.beam_angle,
    #                                 adcpo.config.beam_pattern == 'convex')

    u_inst = \
       beamvel[:,0,:] * rotmat[0, 0] +\
       beamvel[:,1,:] * rotmat[0, 1] +\
       beamvel[:,2,:] * rotmat[0, 2] +\
       beamvel[:,3,:] * rotmat[0, 3]

    v_inst = \
       beamvel[:,0,:] * rotmat[1, 0] +\
       beamvel[:,1,:] * rotmat[1, 1] +\
       beamvel[:,2,:] * rotmat[1, 2] +\
       beamvel[:,3,:] * rotmat[1, 3]

    w_inst = \
       beamvel[:,0,:] * rotmat[2, 0] +\
       beamvel[:,1,:] * rotmat[2, 1] +\
       beamvel[:,2,:] * rotmat[2, 2] +\
       beamvel[:,3,:] * rotmat[2, 3]

    errvel = \
       beamvel[:,0,:] * rotmat[3, 0] +\
       beamvel[:,1,:] * rotmat[3, 1] +\
       beamvel[:,2,:] * rotmat[3, 2] +\
       beamvel[:,3,:] * rotmat[3, 3]

    #ca = rotmat[0,0]
    #b = rotmat[2,0]
    #d = rotmat[3,0]

    #u_inst = ca * (beamvel[:,1,:] - beamvel[:,0,:])
    #v_inst = ca * (beamvel[:,3,:] - beamvel[:,2,:])
    #w_inst = b * beamvel.sum(axis=1)
    #errvel = d *(beamvel[:,1,:] + beamvel[:,0,:] -beamvel[:,3,:] - beamvel[:,2,:]) 

    ## Testing
    #u_inst =beamvel[:,0,:] 
    #v_inst =beamvel[:,1,:] 
    #w_inst =beamvel[:,2,:] 
    #errvel =beamvel[:,3,:] 

    # multiply by the transpose of the rotation matrix
    #u_inst = \
    #   beamvel[:,0,:] * rotmat[0, 0] +\
    #   beamvel[:,1,:] * rotmat[1, 0] +\
    #   beamvel[:,2,:] * rotmat[2, 0] +\
    #   beamvel[:,3,:] * rotmat[3, 0]

    #v_inst = \
    #   beamvel[:,0,:] * rotmat[0, 1] +\
    #   beamvel[:,1,:] * rotmat[1, 1] +\
    #   beamvel[:,2,:] * rotmat[2, 1] +\
    #   beamvel[:,3,:] * rotmat[3, 1]

    #w_inst = \
    #   beamvel[:,0,:] * rotmat[0, 2] +\
    #   beamvel[:,1,:] * rotmat[1, 2] +\
    #   beamvel[:,2,:] * rotmat[2, 2] +\
    #   beamvel[:,3,:] * rotmat[3, 2]

    #errvel = \
    #   beamvel[:,0,:] * rotmat[0, 3] +\
    #   beamvel[:,1,:] * rotmat[1, 3] +\
    #   beamvel[:,2,:] * rotmat[2, 3] +\
    #   beamvel[:,3,:] * rotmat[3, 3]

    return u_inst, v_inst, w_inst, errvel



#def inst2earth(adcpo, fixed_orientation=False):
def inst2earth(u_inst, v_inst, w_inst,\
        heading_deg, pitch_deg, roll_deg,\
        orientation='down',
        heading_offset=None, declination=None,
        fixed_orientation=False):
    """Rotate velocities from the instrument to the earth frame.

    The rotation matrix is taken from the Teledyne RDI
    ADCP Coordinate Transformation manual January 2008
    """
    r = roll_deg * deg2rad
    p = np.arctan(np.tan(pitch_deg * deg2rad) * np.cos(r))
    #p = pitch_deg*deg2rad
    h = heading_deg * deg2rad

    if heading_offset is not None:
        h += heading_offset * deg2rad
    if declination is not None:
        h += declination * deg2rad
    if orientation == 'up':
        r += np.pi

    ch = np.cos(h)
    sh = np.sin(h)
    cr = np.cos(r)
    sr = np.sin(r)
    cp = np.cos(p)
    sp = np.sin(p)

    u = (ch * cr + sh * sp * sr) * u_inst +\
                    (sh * cp) * v_inst +\
                    (ch * sr - sh * sp * cr) * w_inst
    v = (-1 * sh * cr + ch * sp * sr) * u_inst +\
                    (ch * cp) * v_inst +\
                    (-1 * sh * sr - ch * sp * cr) * w_inst
    w = (-cp * sr) * u_inst +\
                    sp * v_inst\
                    + (cp * cr) * w_inst

    ## Multiply by the tranpose of M
    #u = (ch*cr + sh*sp*sr) *u_inst+\
    #    (ch*sp*sr - cr*sh) *v_inst+\
    #    (-cp*sr) *w_inst
    #
    #v = (cp*sh) *u_inst+\
    #    (ch*cp) *v_inst+\
    #    (sp) *w_inst

    #w = (ch*sr - cr*sh*sp) *u_inst+\
    #    (-ch*cr*sp - sh*sr) *v_inst+\
    #    (cp*cr) *w_inst

    return u, v, w

##########
# Depth-bin mapping routines
def binmap(beamvel, distance,\
        pitch_deg, roll_deg,):
    """
    Map the beam data onto the appropriate depth layer
    """

    r = roll_deg * deg2rad
    p = np.arctan(np.tan(pitch_deg * deg2rad) * np.cos(r))
    #p = pitch_deg*deg2rad

    # Calculate the distance of each beam
    dist_1 = distance[..., np.newaxis] * np.cos(p)
    dist_2 = distance[..., np.newaxis] * np.cos(-p)
    dist_3 = distance[..., np.newaxis] * np.cos(r)
    dist_4 = distance[..., np.newaxis] * np.cos(-r)

    nz, nt = dist_1.shape

    beamvelnew = np.zeros_like(beamvel)
    # Build the interpolation objects
    def interp(d1, vel, dist):
        F = interp1d(d1, vel, axis=0, kind='nearest', bounds_error = False)
        return F(dist)

    print('Re-mapping beam data...')
    for ii in range(nt):
        beamvelnew[:,0,ii] = interp(dist_1[:,ii], beamvel[:,0,ii], distance)
        beamvelnew[:,1,ii] = interp(dist_2[:,ii], beamvel[:,1,ii], distance)
        beamvelnew[:,2,ii] = interp(dist_3[:,ii], beamvel[:,2,ii], distance)
        beamvelnew[:,3,ii] = interp(dist_4[:,ii], beamvel[:,3,ii], distance)

    print('Done')
    
    return beamvelnew

