"""
Regrid ROMS model data in 3D
"""

import xarray as xr
import numpy as np
import copy

from sfoda.utils.interpXYZ import interpXYZ
from sfoda.roms.romsio import get_depth
from sfoda.utils.interp1d_numba import interp1d_numba

class ROMSRegridder(object):
    
    grdvar = 'rho'
    interpmethod = 'barycentric'
    romsgrid = None
    P = None
    extrapolate = True
    
    def __init__(self, romsds, x, y, z, **kwargs):
        self.__dict__.update(kwargs)
        
        self.roms = romsds
        assert x.ndim == 1
        assert y.ndim == 1
        

        self.X, self. Y = np.meshgrid(x, y)
        self.x = x
        self.y = y
        self.z = z
        
        self.ny, self.nx = self.X.shape
        self.nz = z.shape[0]
        
        if self.romsgrid is None:
            self.romsgrid = romsds
        
        self.F, self.Ny, self.Nx, self.xyout, self.mask2d, self.maskout = self._build_2d_interp(self.grdvar)
        
        # Get the depth at each point
        if self.grdvar is not 'rho':
            # Depths are stored at rho points
            Fd, _,_,_,maskin,maskd = self._build_2d_interp('rho')
            
        else:
            Fd = self.F
            maskd = self.maskout
            maskin = self.mask2d
        
        self.depth = Fd(self.romsgrid.h.values[maskin]).reshape(self.ny,self.nx)
        self.depth[maskd] = 0

        # U/V rotation angle
        self.uvangle = Fd(self.romsgrid.angle.values[maskin]).reshape(self.ny,self.nx)
        self.uvangle[maskd] = 0


        

        
    def __call__(self, varname, tstep):
        
        # Load the data
        myroms = self.roms.isel(ocean_time=tstep)
        data = myroms[varname].values 

        # Calculate the depth
        zroms = get_depth(myroms.s_rho.values,\
            myroms.Cs_r.values, myroms.hc.values, self.depth,
                          zeta=None, Vtransform=myroms.Vtransform.values)
        
        Nz = zroms.shape[0]
        
        # Interpolate in the horizontal direction
        outdata = np.zeros((Nz, self.ny, self.nx))
        outdataz = np.zeros((self.nz, self.ny, self.nx))

        for kk in range(Nz):
            outdata[kk,...] = self.F(data[kk,self.mask2d]).reshape(self.ny,self.nx)
            # Apply the mask
            outdata[kk,self.maskout] = 0
            
        # Interpolate in the vertical direction
        for ii in range(self.nx):
            for jj in range(self.ny):
                outdataz[:,jj,ii] = interp1d_numba(self.z, zroms[:,jj,ii], outdata[:,jj,ii])
                

        # Loop through all depths and extrapolate top and bottom
        for kk in range(self.nz):
            # bottom cells
            bidx = self.z[kk] <= zroms[0,...]
            if self.extrapolate:
                outdataz[kk, bidx] = outdata[0,bidx]
            else:
                outdataz[kk, bidx] = 0

            # top cells
            tidx = self.z[kk] > zroms[-1,...]
            outdataz[kk, tidx] = outdata[-1,tidx]
        
        return self._to_da(copy.deepcopy(outdataz), varname, tstep)

    def _build_2d_interp(self, grdvar):
        
        # Load the grid
        lonroms = self.romsgrid['lon_{}'.format(grdvar)].values
        latroms = self.romsgrid['lat_{}'.format(grdvar)].values
        mask2d = self.romsgrid['mask_{}'.format(grdvar)].values.astype('bool')
        #print(self.mask2d.shape, lonroms.shape, latroms.shape)
        Ny, Nx = lonroms.shape

        # Reproject
        if self.P is not None:
            Lo, La = self.P.to_xy(lonroms, latroms)
        else:
            Lo, La = lonroms, latroms
        
        # Create the interpolate objects
        xy = np.column_stack((Lo[mask2d],La[mask2d]))
        xyout = np.column_stack([self.X.ravel(), self.Y.ravel()])
        F = interpXYZ(xy, xyout, method=self.interpmethod)
        
        # Output grid mask
        Fmask = interpXYZ(np.column_stack((Lo.ravel(),La.ravel())), xyout, method='nn')
        maskout = Fmask(~mask2d.ravel()).reshape((self.ny,self.nx))
        
        return F, Ny, Nx, xyout, mask2d, maskout
    
    def _to_da(self, data, varname, tstep):
        
        return xr.DataArray(data[None,...], 
             dims=('time','z','y','x',),
             coords={'z':self.z, 'x':self.x, 'y':self.y,
                    'time':np.array([self.roms.ocean_time.values[tstep]])},
                    attrs=self.roms[varname].attrs).chunk({'x':-1,'y':-1,'z':-1})

        
