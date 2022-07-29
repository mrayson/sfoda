"""
Ocean model to SUNTANS boundary/initial condition interpolation code

Usage:

`oceanmodel2bdy(bnd, ds, setUV=True, seth=False, **interp_args)`

or

`oceanmodel2ic(ic, ds, ...)`


"""

from sfoda.utils.interpXYZ import interpXYZ
from sfoda.utils import othertime

from scipy import interpolate
import numpy as np


###
# Main class

class OceanModelInterp(object):
    """
    Class for intperpolating 4D ocean model output stored in a zarr file
    
    Data must be stored on a rectilinear grid with a fixed (z) vertical coordinate
    """
    
    utmzone = 15
    isnorth = True
    
    # Interpolation options
    interpmethod='idw' # 'nn', 'idw', 'kriging', 'griddata'
    NNear=3
    p = 1.0 #  power for inverse distance weighting
    # kriging options
    varmodel = 'spherical'
    nugget = 0.1
    sill = 0.8
    vrange = 250.0

    # Vertical and temporal interpolation methods (scipy interp1d)
    zinterp='linear'
    tinterp='linear'
   
    xdim = 'x'
    ydim = 'y'
    zdim = 'z'
    tdim = 'time'

    rotate_uv = True

    def __init__(self, xr_ds, xi, yi, zi, timei, **kwargs):
        
        self.__dict__.update(kwargs)
        
        self.xi = xi
        self.yi = yi
        self.zi = zi
        self.timei = timei
        
        self.xy_out = np.vstack((xi.ravel(),yi.ravel())).T
        
        self.ds = xr_ds.sel(time=slice(timei[0],timei[-1]))
        
        self.time = self.ds[self.tdim].values
        
        # No mask...
        self.mask_rho = self.ds.mask.values == 0
        
        X,Y = np.meshgrid(self.ds[self.xdim],
                          self.ds[self.ydim])

        self.xy_in = np.vstack((X[self.mask_rho], Y[self.mask_rho])).T
        
        self.Frho = interpXYZ(self.xy_in, self.xy_out,\
            method=self.interpmethod, NNear=self.NNear,\
            p=self.p, varmodel=self.varmodel, nugget=self.nugget,\
                              sill=self.sill,vrange=self.vrange)
        
        
        # Read the vertical coordinate
        self.z = self.ds[self.zdim]
        # Dimesions sizes
        self.Nx = self.xy_out.shape[0]
        self.Nz = self.zi.shape[0]
        self.Nt = len(self.timei)
        
        self.Nz_roms = self.z.shape[0]

        self.Nt_roms = self.time.shape[0]
           
    def interp(self, setUV=True, seth=True):
        """
        Performs the interpolation in this order:
            1) Interpolate onto the horizontal coordinates
            2) Interpolate onto the vertical coordinates
            3) Interpolate onto the time coordinates
        """
        
        zinterp = self.zinterp
        tinterp = self.tinterp
        
        # Initialise the output arrays @ roms time step
        temproms = np.zeros((self.Nt_roms,self.Nz,self.Nx))
        saltroms = np.zeros((self.Nt_roms,self.Nz,self.Nx))
        uroms = np.zeros((self.Nt_roms,self.Nz,self.Nx))
        vroms = np.zeros((self.Nt_roms,self.Nz,self.Nx))

        
        tempold = np.zeros((self.Nz_roms,self.Nx))
        saltold = np.zeros((self.Nz_roms,self.Nx))
        uold = np.zeros((self.Nz_roms,self.Nx))
        vold = np.zeros((self.Nz_roms,self.Nx))

        # Interpolate h
        #h = self.Frho(self.h[self.mask_rho==1])
        
        # Loop through each time step            
        for tstep in range(0,self.Nt_roms):
            print(self.ds[self.tdim].values[tstep],self.ds[self.tdim].values[self.Nt_roms-1],tstep)
            # Read all variables
            self.temp = self.ds['temp'].isel(time=tstep).values
            self.salt = self.ds['salt'].isel(time=tstep).values
            if setUV:
                self.u = self.ds['u'].isel(time=tstep).values
                self.v = self.ds['v'].isel(time=tstep).values
 
                    
            # Interpolate zeta
            #if seth:
            # Always load zeta as it's needed for the depth calculation
            #zetaroms[tstep,:] = self.Frho(self.zeta[self.mask_rho==1])
            
            # Interpolate other 3D variables
            for k in range(0,self.Nz_roms):
                tmp = self.temp[k,...]
                #tempold[k,:] = self.Frho(tmp)
                tempold[k,:] = self.Frho(tmp[self.mask_rho])
                
                tmp = self.salt[k,...]
                saltold[k,:] = self.Frho(tmp[self.mask_rho])
                
                if setUV:
                    tmp = self.u[k,...]
                    uold[k,:] = self.Frho(tmp[self.mask_rho])
                    
                    tmp = self.v[k,...]
                    vold[k,:] = self.Frho(tmp[self.mask_rho])
    
            # Calculate depths (zeta dependent)
            #zroms = get_depth(self.s_rho,self.Cs_r,self.hc, h, zeta=zetaroms[tstep,:], Vtransform=self.Vtransform)
            zroms = self.z

    
            # Interpolate vertically
            for ii in range(0,self.Nx):
                y = tempold[:,ii]
                fillval = (y[0], y[-1]) # Constant end points
                #fillval = y[0]
                #fillval = 'extrapolate'
                Fz = interpolate.interp1d(zroms,y,kind=zinterp,bounds_error=False,fill_value=fillval)
                #Fz = interpolate.PchipInterpolator(zroms[:,ii],y, extrapolate=True)
                temproms[tstep,:,ii] = Fz(self.zi)
                
                y = saltold[:,ii]
                fillval = (y[0], y[-1]) # Constant end points
                Fz = interpolate.interp1d(zroms,y,kind=zinterp,bounds_error=False,fill_value=fillval)
                #Fz = interpolate.PchipInterpolator(zroms[:,ii],y, extrapolate=True)
                saltroms[tstep,:,ii] = Fz(self.zi)
                
                if setUV:
                    y = uold[:,ii]
                    fillval = (y[0], y[-1]) # Constant end points
                    Fz = interpolate.interp1d(zroms,y,kind=zinterp,bounds_error=False,fill_value=fillval)
                    #Fz = interpolate.PchipInterpolator(zroms[:,ii],y, extrapolate=True)
                    uroms[tstep,:,ii] = Fz(self.zi)
                    
                    y = vold[:,ii]
                    fillval = (y[0], y[-1]) # Constant end points
                    Fz = interpolate.interp1d(zroms,y,kind=zinterp,bounds_error=False,fill_value=fillval)
                    #Fz = interpolate.PchipInterpolator(zroms[:,ii],y, extrapolate=True)
                    vroms[tstep,:,ii] = Fz(self.zi)
                    
                
            # End time loop
        
        # Initialise the output arrays @ output time step
        
        # Interpolate temporally
        if self.Nt_roms > 1:
            print('Temporally interpolating ROMS variables...')
            troms = othertime.SecondsSince(self.time)
            tout = othertime.SecondsSince(self.timei)
            #if seth:
            #    print('\tzeta...')
            #    Ft = interpolate.interp1d(troms,zetaroms,axis=0,kind=tinterp,bounds_error=False)
            #    zetaout = Ft(tout)
            #else:
            zetaout=-1

            print('\ttemp...')
            Ft = interpolate.interp1d(troms,temproms,axis=0,kind=tinterp,bounds_error=False)
            tempout = Ft(tout)
            print('\tsalt...')
            Ft = interpolate.interp1d(troms,saltroms,axis=0,kind=tinterp,bounds_error=False)
            saltout = Ft(tout)
            if setUV:
                print('\tu...')
                Ft = interpolate.interp1d(troms,uroms,axis=0,kind=tinterp,bounds_error=False)
                uout = Ft(tout)
                print('\tv...')
                Ft = interpolate.interp1d(troms,vroms,axis=0,kind=tinterp,bounds_error=False)
                vout = Ft(tout)
            else:
                uout = vout = -1
        else:
            zetaout = -1
            tempout = temproms
            saltout = saltroms
            uout = uroms
            vout = vroms
        
        return zetaout, tempout, saltout, uout, vout
    
def oceanmodel2bdy(bnd, ds, setUV=True, seth=False, **interp_args):
    """
    Interpolate to a SUNTANS boundary condition object
    """

    if bnd.N3>0:
        print(72*'#')
        print('Interpolating Type-3 boundary points...')
        F = OceanModelInterp(ds, bnd.xv, bnd.yv, -bnd.z, bnd.time, **interp_args)
        h, T, S, uc, vc = F.interp()

        bnd.T+=T
        bnd.S+=S

        if seth:
            bnd.h+=h
        if setUV:
            bnd.uc+=uc
            bnd.vc+=vc
            
    if bnd.N2>0:
        print(72*'#')
        print('Interpolating Type-2 boundary points...')
        Fe = OceanModelInterp(ds, bnd.xe, bnd.ye, -bnd.z, bnd.time, **interp_args)
        h, T, S, uc, vc = Fe.interp()

        bnd.boundary_T+=T
        bnd.boundary_S+=S

        if setUV:
            bnd.boundary_u += uc
            bnd.boundary_v += vc
            
def oceanmodel2ic(IC, ds, setUV=False, seth=False, **interp_args):
    """
    Interpolate to a SUNTANS boundary condition object
    """

    F = OceanModelInterp(ds, IC.xv, IC.yv, -IC.z_r, np.array([IC.time]), **interp_args)
    _, IC.T, IC.S, IC.uc, IC.vc = F.interp()
    
    if not setUV:
        IC.uc *= 0 
        IC.vc *= 0
    if not seth:
        IC.h *= 0


