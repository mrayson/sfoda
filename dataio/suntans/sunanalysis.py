"""
Tools for in-depth analysis of SUNTANS output

Includes:
    - Volume and tracer budget calculations
    - ...

Matt Rayson
Stanford University
March 2014
"""

from .sunpy import Spatial
from .sunslice import MultiSliceEdge
import soda.utils.mypandas as mpd
from soda.utils.timeseries import timeseries
from soda.utils.maptools import maskShpPoly

import numpy as np
import matplotlib.pyplot as plt
from datetime import timedelta
from netCDF4 import Dataset
from scipy import sparse
import os
import pandas as pd

import pdb

# Global constants
RHO0 = 1000.
Cp = 4186 # specific heat of seawater
GRAV = 9.81

class Energetics(Spatial):

    fluxvar = 'U_F' # U or U_F
    etavar='eta_avg' # 'eta_avg' or 'eta'

    verbose=False
    
    def __init__(self,ncfile,**kwargs):
        """
        Calculates the energy variables from suntans output
        """
        # Initialize the spatial class
        Spatial.__init__(self,ncfile,klayer=[-99])

    def __call__(self,tstep,cellindex=None):
        """
        Calculate the terms for tstep
        """
        if self.verbose: print('Calculating energy at time step: %d'%tstep)
        if cellindex==None: 
            self.cellindex=list(range(self.Nc))
        else:
            self.cellindex=cellindex

        self.tstep=[tstep]

        ###
        # Step 1: Load the flux variable and the vertical depths
        # These are needed for depth integrals and upwind calculations
        ###
        if self.verbose: print('Loading raw model data...')

        self.dzf = self.loadData(variable='dzf')
        # dzf is calculated using max free surface height
        self.dzz = self.loadData(variable='dzz')

        self.u=self.loadData(variable=self.fluxvar)
        if self.fluxvar=='U':
            if self.verbose: print('Calculating U to flow rate...')
            #TBC

        # Load the cell variable used by others at all depth
        self.eta = self.loadData(variable=self.etavar)
        self.uc = self.loadData(variable='uc')
        self.vc = self.loadData(variable='vc')
        self.buoyancy = self.loadData(variable='buoyancy')
        self.nu_v = self.loadData(variable='nu_v')
        if self.hasVar('kappa_tv'):
            self.kappa_tv = self.loadData(variable='kappa_tv')
        else:
            self.kappa_tv = self.nu_v


        # Make sure that all variables = 0 in masked regions...
        # (mask does not work with all operations)
        self.u[self.u.mask]=0
        self.uc[self.uc.mask]=0
        self.vc[self.vc.mask]=0
        self.buoyancy[self.buoyancy.mask]=0
        self.nu_v[self.nu_v.mask]=0
        self.kappa_tv[self.kappa_tv.mask]=0

        # Put all of the terms in a dictionary called... energy
        self.energy={}
        ###
        # Term: Vertical PE flux
        if self.verbose: print('Calculating vertical buoyancy flux...')
        self.energy.update({'B_flux':self.calc_buoyflux()})

        ###
        # Term: Wind work 
        if self.verbose: print('Calculating the wind work...')
        self.energy.update({'W_work':self.calc_windwork()})

        ###
        # Depth integrated KE and PE
        if self.verbose: print('Calculating energy...')
        self.KE = self.calc_KE(u=self.uc,v=self.vc)
        self.energy.update({'KE':self.depthint(self.KE,dz=self.dzz)})

        self.PE = self.calc_PE(b=self.buoyancy)
        self.energy.update({'PE':self.depthint(self.PE,dz=self.dzz)})

        ###
        # Dissipation
        if self.verbose: print('Calculating dissipation...')
        self.energy.update({'diss':self.calc_dissipation()})

        ###
        # Flux terms
        if self.verbose: print('Calculating energy flux divergence terms...')

        # Pressure work flux

        self.energy.update({'uKE':self.calc_fluxdivergence(self.KE)})
        self.energy.update({'uP':self.calc_Pworkflux()})
        self.energy.update({'uPE':self.calc_fluxdivergence(self.PE)})

        # Tide only estimate
        self.energy.update({'ueta':self.calc_fluxdivergence2d(-self.eta*GRAV)})


    def write2netcdf(self,outfile,trange):
        """
        Write all time steps in trange to an output file

        !! Note that all terms are converted to Wm-2 (multiplied by rho0) !!
        !! Divergent terms are divided by cell area (self.Ac) !!!
        """
        tstep = list(range(0,self.Nt))[trange[0]:trange[1]]
        # Write the output to netcdf
        print('Writing the output to netcdf...')

        self.writeNC(outfile)

        nc = Dataset(outfile,'a')
        nc.Title = 'SUNTANS energy output'


        nc.close()

        # Create the new variable names
        self.create_nc_var(outfile, 'time', ('time',),\
            {'long_name':'time','units':'seconds since 1990-01-01 00:00:00'})

        self.create_nc_var(outfile, 'KEz', ('time','Nc'),\
            {'long_name':'Depth-integrated kinetic energy',\
            'units':'J m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'PEz', ('time','Nc'),\
            {'long_name':'Depth-integrated potential energy',\
            'units':'J m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'uP', ('time','Nc'),\
            {'long_name':'Depth-integrated pressure work divergence',\
            'units':'W m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'uKE', ('time','Nc'),\
            {'long_name':'Depth-integrated kinetic energy flux divergence',\
            'units':'W m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'uPE', ('time','Nc'),\
            {'long_name':'Depth-integrated potential energy flux divergence',\
            'units':'W m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'ueta', ('time','Nc'),\
            {'long_name':'Depth-integrated tidal energy flux divergence',\
            'units':'W m-2','coordinates':'yv xv'})

        self.create_nc_var(outfile, 'W_work', ('time','Nc'),\
            {'long_name':'Wind work',\
            'units':'W m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'B_flux', ('time','Nc'),\
            {'long_name':'Turbulent vertical buoyancy flux (KE->PE)',\
            'units':'W m-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'diss', ('time','Nc'),\
            {'long_name':'Dissipation rate',\
            'units':'W m-2','coordinates':'yv xv'})
        
        # Testing variables
        self.create_nc_var(outfile, 'S2', ('time','Nk','Nc'),\
            {'long_name':'Shear squared',\
            'units':'s-2','coordinates':'yv xv'})
        self.create_nc_var(outfile, 'Pressure', ('time','Nk','Nc'),\
            {'long_name':'Pressure',\
            'units':'Pa','coordinates':'yv xv'})
 

 
        
        # Calculate the energy for each time step and write the output
        print('Writing the variable data to netcdf...')
        nc = Dataset(outfile,'a')

        for ii, tt in enumerate(tstep):
            # Call the object to calculate the variables
            print('Writing energy for timestep %d of %d...'%(tt,tstep[-1]))
            self.__call__(tt)

            # Write the variable data out
            nc.variables['time'][ii]=self.timeraw[tt]

            nc.variables['KEz'][ii,:]=self.energy['KE']*RHO0
            nc.variables['PEz'][ii,:]=self.energy['PE']*RHO0
            nc.variables['uP'][ii,:]=self.energy['uP']/self.Ac*RHO0
            nc.variables['uKE'][ii,:]=self.energy['uKE']/self.Ac*RHO0
            nc.variables['uPE'][ii,:]=self.energy['uPE']/self.Ac*RHO0
            nc.variables['ueta'][ii,:]=self.energy['ueta']/self.Ac*RHO0
            nc.variables['W_work'][ii,:]=self.energy['W_work']*RHO0
            nc.variables['B_flux'][ii,:]=self.energy['B_flux']*RHO0
            nc.variables['diss'][ii,:]=self.energy['diss']*RHO0
            # Testing variables
            nc.variables['S2'][ii,:,:]=self.S2
            nc.variables['Pressure'][ii,:,:]=self.pressure*RHO0

        nc.close()

    def gradZ(self,phi,dzz,dzmin=0.01):
        """
        Overloaded vertical gradient calculation function

        Make sure the calculation is consistent with turbulence.c 
        Gradients are evaluated at k-1/2

        """
        Nc = phi.shape[1]
        dphi_dz = np.zeros((self.Nkmax+1,Nc))

        #dzz values less than dzmin are set to dzmin
        dzz[dzz<dzmin]=dzmin

        # Calculate mid-point gradients
        dphi_dz[1:-1,:] = 2.0 * (phi[0:-1,:] - phi[1:,:])/ \
            (dzz[0:-1,:]+dzz[1:,:])

        # Specify the surface gradient the same as the next layer
        ctop = self.getctop(self.eta)
        j = list(range(Nc))
        dphi_dz[ctop[j],j] = dphi_dz[ctop[j]+1,j]

        # Specify the seabed gradients
        dphi_dz[self.Nk[j]+1,j]=dphi_dz[self.Nk[j],j]

        # Return the average at the mid-layer depth
        return 0.5*(dphi_dz[1:,:] + dphi_dz[0:-1,:])

    def calc_fluxdivergence(self,phi):
        """
        Calculates the flux divergece of a cell-centered scalar, phi.
        """
        # Map the data onto the edge
        phi_e = np.zeros((self.Nkmax,self.Ne))
        for k in range(self.Nkmax):
            phi_e[k,:] = \
                self.get_edgevar(phi[k,:],k=k,U=self.u[k,:],method='upwind')

        face = self.face.copy()
        normal = 1.0*self.normal # Convert to float
        mask = face.mask.copy()

        # Create a mask so that tthe masked face values are not included in the
        # flucx calculation
        facemask = np.zeros((self.Nc,self.maxfaces))
        facemask[mask==False]=1.0

        face[mask]=0 # Index the mask to the first cell 
        #   (this is multiplied by zero later..)

        # Calculate the fluxes at all cells - dimensions: [Nk, Nc, nfaces]
        flux_cell = phi_e[...,face] * self.u[...,face] * normal * facemask


        # Sum along all faces - dimensions: [Nk, Nc]
        flux_div = flux_cell.sum(axis=-1)


        # Return the depth integrated divergence
        return flux_div.sum(axis=0)

    def calc_fluxdivergence2d(self,phi):
        """
        Calculates the flux divergece of a cell-centered scalar, phi.
        """
        #depth-integrate the flux rate
        U = np.sum(self.u,axis=0)
        de = self.get_edgevar(self.dv,method='max')
        U/=de # Divide by the edge depth (u_bar)

        # Map the data onto the edge
        phi_e = self.get_edgevar(phi,k=0,U=U,method='upwind')

        face = self.face.copy()
        normal = 1.0*self.normal # Convert to float
        mask = face.mask.copy()

        # Create a mask so that tthe masked face values are not included in the
        # flucx calculation
        facemask = np.zeros((self.Nc,self.maxfaces))
        facemask[mask==False]=1.0

        face[mask]=0 # Index the mask to the first cell 
        #   (this is multiplied by zero later..)

        # Calculate the fluxes at all cells - dimensions: [Nc, nfaces]
        flux_cell = phi_e[face] * U[face] * normal * facemask


        # Sum along all faces - dimensions: [Nc]
        return flux_cell.sum(axis=-1)


    def calc_Pworkflux(self):
        """
        Calculate the pressure work flux divergence for all grid cells
        """

        # Calculate pressure at the mid-point
        # Note that this is already normalized by rho0
        #rho = self.buoyancy/GRAV*RHO0+RHO0
        #self.pressure = self.depthint(-GRAV*rho,dz=self.dzz,cumulative=True)
        # Buoyancy only
        self.pressure = self.depthint(self.buoyancy,dz=self.dzz,cumulative=True)

        # Need to add the free-surface contribution???
        # Shouldn't be necessary since dzz varies with the free-surface
        #H = self.depthint(self.dzz,dz=self.dzz,cumulative=True) # total depth
        #self.pressure += H*GRAV # H = eta - z
        self.pressure += self.eta*GRAV
        
        #return self.calc_fluxdivergence(self.pressure/RHO0)
        return self.calc_fluxdivergence(self.pressure)

    def calc_windwork(self):
        """
        Calculate the wind work component
        """
        u_surf = self.get_surfacevar(self.uc,self.eta)
        v_surf = self.get_surfacevar(self.vc,self.eta)
        tau_x = self.loadData(variable='tau_x')
        tau_y = self.loadData(variable='tau_y')

        return (u_surf*tau_x + v_surf*tau_y)/RHO0
        
    def calc_buoyflux(self):
        """
        Calculates the vertical flux of buoyancy:

            B_f = K_v * db/dz

        Returns the depth-integral.
        """

        db_dz = self.gradZ(self.buoyancy,self.dzz)

        return self.depthint(self.kappa_tv*db_dz,dz=self.dzz)

    def calc_dissipation(self):
        r"""
        Calculates the depth-integrated dissipation
            eps = nu_v * (du/dz^2 + dv/dz^2)
        """
        du_dz = self.gradZ(self.uc,self.dzz)
        dv_dz = self.gradZ(self.vc,self.dzz)
        self.S2 = (du_dz**2 + dv_dz**2)
        # Zero the seabed shear - it is too large??
        self.S2[self.Nk,list(range(self.Nc))]=0

        diss = self.nu_v * self.S2

        return self.depthint(diss,dz=self.dzz) 

########################
########################
def energy_budget(energyfile,polyfile,trange):
    """
    # Area-integrate the energy terms
    """
    varnames = ['KEz','PEz','uP','uKE','uPE','ueta','W_work','B_flux','diss']

    # Load the energy file as a suntans object
    sun = Spatial(energyfile)

    # Create the mask
    mask,maskpoly = maskShpPoly(sun.xv,sun.yv,polyfile)

    # Initialise the output dictionary
    tstep = list(range(0,sun.Nt))[trange[0]:trange[1]]
    nt = len(tstep)

    budget ={}
    for vv in varnames:
        budget.update({vv:np.zeros((nt,))})

    for ii,tt in enumerate(tstep):
        print('Area-integrating step: %d of %d...'%(ii,tstep[-1]))
        for vv in varnames:
            sun.tstep=[tt]
            data = sun.loadData(variable=vv)
            budget[vv][ii],areatotal = sun.areaint(data,mask)

    budget.update({'time':sun.time[tstep]})

    # Calculate the time-rate of change of KE and PE
    dt = sun.timeraw[1]-sun.timeraw[0]
    budget.update({'dKE_dt':np.zeros((nt,))})
    budget.update({'dPE_dt':np.zeros((nt,))})
    budget['dKE_dt'][1::] = (budget['KEz'][1::]-budget['KEz'][0:-1])/dt
    budget['dPE_dt'][1::] = (budget['PEz'][1::]-budget['PEz'][0:-1])/dt

    return budget


########################
########################
def calc_avg_budget(sun, trange, cellindex,plot=False):
    """
    Calculate the volume, temperature and salt budgets from
    an average output file. 

    These calculations are very specific to the variables
    stored in the averages file.
    """

    # Load the SUNTANS average file object
    sun.klayer=[-99]
    #sun = Spatial(avgfile,klayer=[-99])

    # Calculate the time dimensions
    tstep = list(range(0,sun.Nt))[trange[0]:trange[1]]
    nt = len(tstep)
    time = sun.time[tstep]

    dt = sun.globalatts['dt']*sun.globalatts['ntaverage']

    # Remove cells that are next to type-2 or 3 edges here
    # ...
    #facemark=sun.get_facemark()
    #for cc in cellindex:
    #    if facemark[cc] in [2,3]:
    #        print 'Removing edge cell index = %d'%cc
    #        cellindex.remove(cc)

    Nc = len(cellindex)

    # Calculate some grid variables
    area = sun.Ac[cellindex]
    sumarea = np.sum(area)

    face = sun.face[cellindex,:] # edge pointers for each cell
    normal = 1.0*sun.normal[cellindex,:] 

    # Create a mask so that the masked face values are not included
    # in the flux calculations
    facemask = np.zeros_like(normal)
    facemask[face.mask==False]=1.0

    face[face.mask]=0 # Index masked cells to the first cell (this is multiplied
        # by zero

    # Initialise the output variables
    # Sum of fluxes
    Mass_f = np.zeros((nt,),np.float)
    Salt_f = np.zeros((nt,),np.float)
    Temp_f = np.zeros((nt,),np.float)

    # Volume integrals
    V = np.zeros((nt,),np.float)
    s_V = np.zeros((nt,),np.float)
    T_V = np.zeros((nt,),np.float)

    # Surface fluxes (T and S only)
    s_surf = np.zeros((nt,),np.float)
    T_surf = np.zeros((nt,),np.float)


    ###
    # Start stepping through and read all variable time step by time step
    ###
    for ii,tt in enumerate(tstep):
        sun.tstep=[tt]
        print('Calculating budget for time = %d of %d'%(tt,tstep[-1]))
        # Load the depth-average and flux quantities
        s_dz = sun.loadDataRaw(variable='s_dz')
        T_dz = sun.loadDataRaw(variable='T_dz')
        eta = sun.loadDataRaw(variable='eta')

        Sflux = sun.loadDataRaw(variable='s_F')
        Tflux = sun.loadDataRaw(variable='T_F')
        Mflux = sun.loadDataRaw(variable='U_F') #[m3 s-1] * [s] = [m3]

        # Subset the variables at cell index only
        eta = eta[cellindex]
        s_dz = s_dz[cellindex]
        T_dz = T_dz[cellindex]

        # Calculate the fluxes for each cell [Nk, Nc, maxfaces]
        Mflux_cell = Mflux[...,face] * normal * facemask
        Sflux_cell = Sflux[...,face] * normal * facemask
        Tflux_cell = Tflux[...,face] * normal * facemask

        # Compute the total mass/tracer flux in/out of each cell
        # sum along all dimension edges
        Mass = Mflux_cell.sum(axis=-1)  
        Salt = Sflux_cell.sum(axis=-1)  
        Temp = Tflux_cell.sum(axis=-1)
        # Sum along all depth
        Mass = Mass.sum(axis=0)
        Salt = Salt.sum(axis=0)
        Temp = Temp.sum(axis=0)

        # Sum all cells 
        Mass_f[ii] = Mass.sum()
        Salt_f[ii] = Salt.sum()
        Temp_f[ii] = Temp.sum()

        # Calculate the volume integrals
        s_V[ii] = np.sum(s_dz*area,axis=-1).squeeze() # m3 S
        T_V[ii] = np.sum(T_dz*area,axis=-1).squeeze() # m3 C
        V[ii] = np.sum(eta*area,axis=-1).squeeze() # m3 [volume]

        # Get the surface temp and salinity  flux arrays
        if sun.hasVar('Hs'):
            # Load the surface flux quantities
            Hs = sun.loadDataRaw(variable='Hs')
            Hsw = sun.loadDataRaw(variable='Hsw')
            Hl = sun.loadDataRaw(variable='Hl')
            Hlw = sun.loadDataRaw(variable='Hlw')

            # Convert heat flux [W m-2] -> temperature flux
            Qs = (Hs+Hl+Hlw+Hsw)/(RHO0*Cp) # units [C m s-1]

            # Surface flux contribution
            T_surf[ii] = np.sum(Qs[...,cellindex]*area) # units [C m3 s-1]
        else:
            T_surf[ii] = 0

        if sun.hasVar('EP'):
            EPs0 = sun.loadDataRaw(variable='EP')
            s_surf[ii] = np.sum(EPs0[...,cellindex]*area)  # units [psu m3 s-1]
        else:
            s_surf[ii] = 0


    ##########
    # units are:
    ##########
    # s_V [ppt m3]
    # T_V [C m3]
    # eta [m3]

    # Mass_f [m3 s-1]
    # Salt_f [ppt m3 s-1]
    # Temp_f [C m3 s-1]

    ###
    # Compute each of the terms in the budget

    # Tendency
    Tend_V = (V[:-1]-V[1:]).squeeze()/dt # m3 s-1
    Tend_s = (s_V[:-1]-s_V[1:]).squeeze()/dt # psu m3 s-1
    Tend_T = (T_V[:-1]-T_V[1:]).squeeze()/dt # C m3 s-1

    # Advective fluxes
    Adv_V = Mass_f[1:]# m3 s-1
    Adv_s = Salt_f[1:]# psu s-1
    Adv_T = Temp_f[1:]# C s-1

    # Surface fluxes (note change of sign)
    Sflux_T = -T_surf[1:]# C m3 s-1
    Sflux_s = s_surf[1:]# psu m3 s-1

    # Compute the error (residual) in each budget
    Err_V =(Tend_V - Adv_V) 
    Err_T = (Tend_T - Adv_T - Sflux_T) 
    Err_s = (Tend_s - Adv_s - Sflux_s)

    # Output time
    time = time[1:]

    # Save the output as a dictionary
    budget = {'time':time,\
              'cellindex':cellindex,\
              'Volume':{'Advection':Adv_V,'Tendency':Tend_V,'Residual':Err_V},\
              'Temp':{'Advection':Adv_T,'Tendency':Tend_T,'Surface_Flux':Sflux_T,'Residual':Err_T},\
              'Salt':{'Advection':Adv_s,'Tendency':Tend_s,'Surface_Flux':Sflux_s,'Residual':Err_s},\
              }

    if plot:
        # Free-surface
        fig1=plt.figure()
        f1ax1=fig1.add_subplot(2,1,1)
        plt.title('Volume budget')
        plt.plot(time,Tend_V,'b',linewidth=2)
        plt.plot(time,Adv_V,'r')
        plt.ylabel('$m^3 \ s^{-1}$')
        plt.ylim(Tend_V.min(),Tend_V.max())
        plt.legend(('Tendency','Advection'))

        ax2=fig1.add_subplot(2,1,2,sharex=f1ax1)
        plt.plot(time,Err_V)
        plt.ylabel('error')

        fig2=plt.figure()
        f2ax1=fig2.add_subplot(2,1,1)
        plt.title('Temperature budget')
        plt.plot(time,Tend_T,'b',linewidth=2)
        plt.plot(time,Adv_T,'r')
        plt.plot(time,Adv_T + Sflux_T,'k')
        plt.grid(b=True)
        plt.ylabel(r'$^\circ C \ m^3 \ s^{-1}$')
        plt.legend(('Tendency','Advection','Adv. + Sflux'))

        f2ax1=fig2.add_subplot(2,1,2,sharex=f2ax1)
        plt.title('Temperature budget')
        plt.plot(time,Err_T)
        plt.ylabel('error')

        fig3=plt.figure()
        f3ax1=fig3.add_subplot(2,1,1)
        plt.title('Salt budget')
        plt.plot(time,Tend_s,'b',linewidth=2)
        plt.plot(time,Adv_s,'r')
        plt.plot(time,Adv_s + Sflux_s,'k')
        plt.grid(b=True)
        plt.ylabel('$psu \ m^3 \ s^{-1}$')
        plt.legend(('Tendency','Advection','Adv. + Sflux'))

        f2ax1=fig3.add_subplot(2,1,2,sharex=f3ax1)
        plt.plot(time,Err_s)
        plt.ylabel('error')

        plt.figure()
        sun.plotmesh()
        plt.plot(sun.xv[cellindex],sun.yv[cellindex],'m.')
        plt.show()

    return budget
    #

def calc_isopycnal_discharge(ncfile,xpt,ypt,saltbins,tstart,tend,scalarvar='salt'):
    """
    Calculates the discharge as a function of salinity along
    a transect, defined by xpt/ypt, in the suntans model

    Returns a dictionary with the relevant variables
    """
    nbins = saltbins.shape[0]

    # Load the slice object and extract the data
    SE = MultiSliceEdge(ncfile,xpt=xpt,ypt=ypt)
    #    if SE==None:
    #        SE = SliceEdge(ncfile,xpt=xpt,ypt=ypt)
    #        SE.tstep = range(SE.Nt)
    #    else:
    #        SE.update_xy(xpt,ypt)
    #
 
    SE.tstep = SE.getTstep(tstart,tend)

    print('Loading the salt flux data...')
    #s_F_all= SE.loadData(variable='s_F')
    s_F_all= SE.loadData(variable=scalarvar)
    print('Loading the flux data...')
    Q_all = SE.loadData(variable='U_F')


    def Q_S_flux(salt,Q,saltbins,normal):
        # mask sure masked values are zeroed
        #s_F[s_F.mask]=0
        #Q[Q.mask]=0
        Q = Q*normal

        Nt,Nk,Ne = Q.shape

        #salt = np.abs(s_F)/np.abs(Q)
        #salt[np.isnan(salt)]=0

        Ns = saltbins.shape[0]
        ds = np.diff(saltbins).mean()

        ###
        # Calculate Q(s,x)
        ###
        # Create an arrayo
        #Nt = len(SE.tstep)
        #ne = len(SE.j) # number of edges
        jindex = np.arange(0,Ne)
        jindex = np.repeat(jindex[np.newaxis,np.newaxis,:],Nt,axis=0)
        jindex = np.repeat(jindex,SE.Nkmax,axis=1)

        # Groups the salt matrix into bins
        sindex = np.searchsorted(saltbins,salt)
        sindex[sindex>=Ns]=Ns-1

        #tindex = np.arange(0,Nt)
        #tindex = np.repeat(tindex[:,np.newaxis,np.newaxis],ne,axis=-1)
        #tindex = np.repeat(tindex,SE.Nkmax,axis=1)

        # Calculate the salt flux for each time step

        Qs = np.zeros((Nt,Ns,Ne))# 
        #Fs = np.zeros((Nt,Ne))# 
        #dQds = np.zeros((Nt,Ns,Ne))# put salt in last dimension for easy  multiplication
        for tt in range(Nt):
            
            # Create an array output array Q_S 
            # This sums duplicate elements
            Q_S = sparse.coo_matrix((Q[tt,...].ravel(),\
                (sindex[tt,...].ravel(),jindex[tt,...].ravel())),\
                shape=(Ns,Ne)).todense()
            Qs[tt,...] = np.array(Q_S)#/Nt # Units m^3/s

            ####
            ## THIS IS WRONG DON'T USE
            ####
            ## Compute the gradient (this gives the same result after
            ## integration)
            #dQ_ds, dQ_de = np.gradient(Qs[tt,...],ds,1)

            ##dQtmp = -1*np.array(dQ_ds).T
            #dQds[tt,...] = -1*dQ_ds
            #
            #Fs[tt,:] = np.sum(-1*dQds[tt,...].T*saltbins ,axis=-1)

        output = {'time':SE.time[SE.tstep],'saltbins':saltbins,\
            'Qs':Qs}
            #'dQds':dQds,'Qs':Qs,'Fs':Fs}

        return output




    def Q_S_flux_old(s_F,Q,saltbins,x,normal,area,dz):
        # mask sure masked values are zeroed
        #s_F[s_F.mask]=0
        #Q[Q.mask]=0
        Q = Q*normal

        Nt,Nk,ne = Q.shape

        salt = np.abs(s_F)/np.abs(Q)
        salt[np.isnan(salt)]=0

        # Calculate the mean Q
        Qbar = np.sum( np.sum(Q,axis=-1),axis=0)/Nt

        ###
        # Calculate Q(s,x)
        ###
        # Create an arrayo
        #Nt = len(SE.tstep)
        #ne = len(SE.j) # number of edges
        jindex = np.arange(0,ne)
        jindex = np.repeat(jindex[np.newaxis,np.newaxis,:],Nt,axis=0)
        jindex = np.repeat(jindex,SE.Nkmax,axis=1)

        # Groups the salt matrix into bins
        sindex = np.searchsorted(saltbins,salt)
        sindex[sindex>=nbins]=nbins-1

        # Create an array output array Q_S 
        # This sums duplicate elements
        Q_S_x = sparse.coo_matrix((Q.ravel(),(sindex.ravel(),jindex.ravel())),\
            shape=(nbins,ne)).todense()
        Q_S_x = np.array(Q_S_x)#/Nt # Units m^3/s

        ###
        # Calculate Q(s,t)
        ###
        # Create an arrayo
        tindex = np.arange(0,Nt)
        tindex = np.repeat(tindex[:,np.newaxis,np.newaxis],ne,axis=-1)
        tindex = np.repeat(tindex,SE.Nkmax,axis=1)


        # Create an array output array Q_S 
        # This sums duplicate elements
        Q_S_t = sparse.coo_matrix((Q.ravel(),(sindex.ravel(),tindex.ravel())),\
            shape=(nbins,Nt)).todense()
        Q_S_t = np.array(Q_S_t)#/ne # Units m^3/s


        ###
        # Calculate Q(s)
        ###
        Q_S = np.bincount(sindex.ravel(),weights=Q.ravel(),minlength=nbins)

        ###
        # Calculate the gradients with respect to S
        ###
        ds = np.diff(saltbins).mean()
        dsdt_inv = 1./(ds*Nt)
        saltbins = 0.5*(saltbins[1:] + saltbins[0:-1]) 

        # Units are: [m^3 s^-1 psu^-1]
        dQ_S_x = np.diff(Q_S_x,axis=0) * dsdt_inv
        dQ_S_t = np.diff(Q_S_t,axis=0) * dsdt_inv
        dQ_S = np.diff(Q_S,axis=0) * dsdt_inv

        ###
        # Now integrate to calculate the flux terms
        # See Macready 2011 eq. 3 and 4
        ###
        ind_in = dQ_S>=0
        ind_out = dQ_S<0
        Fin = np.sum(saltbins[ind_in] * -dQ_S[ind_in]*ds) 
        Fout = np.sum(saltbins[ind_out] * -dQ_S[ind_out]*ds) 
        Qin = np.sum(-dQ_S[ind_in]*ds) 
        Qout = np.sum(-dQ_S[ind_out]*ds) 

        # Put all of the relevant variables into a dictionary
        output = {'x':x,'time':SE.time[SE.tstep],'saltbins':saltbins,\
            'dQ_S':dQ_S,'dQ_S_x':dQ_S_x,'dQ_S_t':dQ_S_t,\
            'F_in':Fin,'F_out':Fout,'Q_in':Qin,'Q_out':Qout}

        return output,Qbar

    output =[]
    outputf =[]
    print('Calculating slice fluxes...')
    ii=-1
    Qbar = np.zeros((len(Q_all),SE.Nkmax))
    for s_F,Q in zip(s_F_all,Q_all):
        ii+=1
        x = SE.slices[ii]['distslice'][1:]
        normal = SE.slices[ii]['normal']
        area = SE.slices[ii]['area']
        dx = SE.slices[ii]['dx']
        #tmp,Qbar[ii,:] = Q_S_flux_old(s_F,Q,saltbins,x,normal,area,SE.dz)
        tmp = Q_S_flux(s_F,Q,saltbins,normal)
        output.append(tmp)

        ## Also calculate the filtered
        TS = timeseries(SE.time[SE.tstep],s_F)
        TS.godinfilt()
        #s_F_filt = TS.y.reshape(s_F.shape)
        s_F_filt = TS.y.T
        TS = timeseries(SE.time[SE.tstep],Q)
        TS.godinfilt()
        #Q_filt = TS.y.reshape(Q.shape)
        Q_filt = TS.y.T
        tmp = Q_S_flux(s_F_filt,Q_filt,saltbins,normal)

        #TS = pd.Panel(s_F, items=SE.time[SE.tstep])
        #TSf = mpd.godin(TS)
        #pdb.set_trace()

        outputf.append(tmp)


    # Calculate the eulerian mean
    #Sbar = SE.mean(s_F_all,axis='depth')
    #Qbar = SE.mean(Q_all,axis='depth')
    #z = SE.z_r

    return output,outputf, SE
    #return output,SE,(Sbar,Qbar,z)

def calc_salt_flux(ncfile,xpt,ypt,tstart,tend):
    return calc_scalar_flux(ncfile,xpt,ypt,tstart,tend)

def calc_scalar_flux(ncfile,xpt,ypt,tstart,tend, \
        scalarvar='salt', fluxvar='s_F', edgemethod=0):
    """
    Calculate the scalar and discharge flux through the sections in xpt,ypt
    """
    SE = MultiSliceEdge(ncfile,xpt=xpt,ypt=ypt, edgemethod=edgemethod)

    #SE = SliceEdge(ncfile,xpt=xpt,ypt=ypt)

    if tstart==None:
        SE.tstep = list(range(len(SE.time)))
    else:
        SE.tstep = SE.getTstep(tstart,tend)

    print('Loading the area data...')
    A_all= SE.loadData(variable='area')
    print('Loading the salt data...')
    S_all= SE.loadData(variable=scalarvar) # psu 
    print('Loading the salt flux data...')
    F_all= SE.loadData(variable=fluxvar) # psu m3 s-1
    print('Loading the flux data...')
    Q_all = SE.loadData(variable='U_F')


    def calc_flux_light(Q_all,A_all,normal):

        Nt,Nk,Ne = Q_all.shape

        # Point Q in the right direction
        Q_all = Q_all * normal

        Q = Q_all.sum(axis=-1).sum(axis=-1)
        A = A_all.sum(axis=-1).sum(axis=-1)

        # Decompose into inflow and outflow
        ind_in = Q_all > 0
        ind_out = ind_in == False

        Qin = np.zeros((Nt,))
        Qout = np.zeros((Nt,))

        for ii in range(Nt):
            Qin[ii] = np.sum(Q_all[ii,ind_in[ii,...]])
            Qout[ii] = np.sum(Q_all[ii,ind_out[ii,...]])

        return Q,Qin,Qout,A


    def calc_flux(Q_all,F_all, S_all,A_all,normal):

        Nt,Nk,Ne = Q_all.shape

        # Point Q in the right direction
        Q_all = Q_all * normal
        F_all = F_all * normal

        Q = Q_all.sum(axis=-1).sum(axis=-1)
        F = F_all.sum(axis=-1).sum(axis=-1)
        A = A_all.sum(axis=-1).sum(axis=-1)
        S = (S_all*A_all).sum(axis=-1).sum(axis=-1)/A

        # Decompose into inflow and outflow
        ind_in = Q_all > 0
        ind_out = ind_in == False

        Qin = np.zeros((Nt,))
        Qout = np.zeros((Nt,))
        Fin = np.zeros((Nt,))
        Fout = np.zeros((Nt,))

        for ii in range(Nt):
            Qin[ii] = np.sum(Q_all[ii,ind_in[ii,...]])
            Qout[ii] = np.sum(Q_all[ii,ind_out[ii,...]])
            Fin[ii] = np.sum(F_all[ii,ind_in[ii,...]])
            Fout[ii] = np.sum(F_all[ii,ind_out[ii,...]])

        return Q,Qin,Qout,F,Fin,Fout,S,A

    data = []
    for ii in range(len(Q_all)):
        Q,Qin,Qout,F,Fin,Fout,S,A =\
            calc_flux(Q_all[ii],F_all[ii],S_all[ii],A_all[ii],SE.slices[ii]['normal'])
        #Q, Qin, Qout, A =\
        #    calc_flux_light(Q_all[ii], A_all[ii], SE.slices[ii]['normal'])

        data.append({'time':SE.time[SE.tstep],\
                'Q':Q,'Qin':Qin,'Qout':Qout,\
                'F':F,'Fin':Fin,'Fout':Fout,\
                'S':S,\
                'area':A})

    return data, SE

def volume_integrate(ncfile,varnames,shpfiles,constantdzz=False):
    """
    Volume integrate a suntans variable for all time in the domain 
    specified with a shpfile polygon
    """
    # Use numexpr to try and speed things up
    import numexpr as ne

    # Load the spatial object
    sun = Spatial(ncfile,klayer=[-99])
    Ac = sun.Ac

    # Calculate the mask region
    if type(shpfiles) != type([]):
        shpfiles = [shpfiles]

    masks = []
    polynames = []
    for shpfile in shpfiles:
        mask, maskpoly = maskShpPoly(sun.xv,sun.yv,shpfile)
        masks.append(mask)
        polynames.append(os.path.splitext(os.path.basename(shpfile))[0])

    # Create a dictionary with output info
    data={}
    for poly  in polynames:
        data.update({poly:{'V':np.zeros((sun.Nt,)),'time':sun.time}})
        for varname in varnames:
            data[poly].update({varname:np.zeros((sun.Nt,))})

    # Fix dzz for testing
    sun.tstep = [0]
    dzz = sun.loadData(variable='dzz')
    h = ne.evaluate("sum(dzz,axis=0)")

    #dzz = np.repeat(sun.dz[:,np.newaxis],sun.Nc,axis=1)
    for ii in range(sun.Nt):
        sun.tstep = [ii]
        print('Volume integrating for time step: %d of %d...'%(ii,sun.Nt))

        # Load the depth and mean age arrays
        if not constantdzz:
            dzz = sun.loadData(variable='dzz')
            # Calculate the total volume 
            #h = np.sum(dzz,axis=0)
            h = ne.evaluate("sum(dzz,axis=0)")
        
        for varname in varnames:
            tmp = sun.loadData(variable=varname)

            for mask,poly in zip(masks,polynames):
                V, A = sun.areaint(h,mask=mask) # Depth*area
                data[poly]['V'][ii] = V

                # Get the depth-integral 
                #tmp_dz = sun.depthint(tmp,dz=dzz)
                tmp_dz = ne.evaluate("sum(tmp*dzz,axis=0)")
                
                # Calculate the volume-integral
                #tmp_dV, A = sun.areaint(tmp_dz,mask=mask)
                tmp_dV = ne.evaluate("sum(tmp_dz*Ac*mask)")

                data[poly][varname][ii]=tmp_dV/V

    return data

def area_integrate(ncfile,varnames,shpfiles):
    """
    Area integrate a suntans variable for all time in the domain 
    specified with a shpfile polygon
    """
    # Use numexpr to try and speed things up
    import numexpr as ne

    # Load the spatial object
    sun = Spatial(ncfile,klayer=[-99])
    Ac = sun.Ac

    # Calculate the mask region
    if type(shpfiles) != type([]):
        shpfiles = [shpfiles]

    masks = []
    polynames = []
    for shpfile in shpfiles:
        mask, maskpoly = maskShpPoly(sun.xv,sun.yv,shpfile)
        masks.append(mask)
        polynames.append(os.path.splitext(os.path.basename(shpfile))[0])

    # Create a dictionary with output info
    data={}
    for poly  in polynames:
        data.update({poly:{'V':np.zeros((sun.Nt,)),'time':sun.time}})
        for varname in varnames:
            data[poly].update({varname:np.zeros((sun.Nt,))})

    sun.tstep = list(range(sun.Nt))
    for varname in varnames:
        print('Area integrating varibles: %s ...'%(varname))
        tmp = sun.loadData(variable=varname)
        for mask,poly in zip(masks,polynames):
            tmp_dA = ne.evaluate("sum(tmp*Ac*mask,axis=1)")
            data[poly][varname]=tmp_dA


    return data

def river_discharge(ncfile,shpfiles):
    """
    Calculates the river flux of all type-2 boundaries located with each polygon 
    specified with a shpfile polygon
    """

    # Load the spatial object
    sun = Spatial(ncfile,klayer=[-99])

    # Indentify the river cells
    ind = np.argwhere(sun.mark==2).ravel()
    #sun.j = ind

    # Calculate the mask region
    if type(shpfiles) != type([]):
        shpfiles = [shpfiles]

    masks = []
    polynames = []
    for shpfile in shpfiles:
        mask, maskpoly = maskShpPoly(sun.xe[ind],sun.ye[ind],shpfile)
        masks.append(mask)
        polynames.append(os.path.splitext(os.path.basename(shpfile))[0])

    # Create a dictionary with output info
    data={}
    for poly  in polynames:
        data.update({poly:{'Q_r':np.zeros((sun.Nt,)),'time':sun.time}})

    print('Loading the data...')
    #U = sun.loadData(variable='U_F')
    sun.tstep = list(range(sun.Nt))
    U = np.zeros((sun.Nt,sun.Nkmax,ind.shape[0]))
    for tt in range(sun.Nt):
        sun.tstep = tt 
        tmp = sun.loadData(variable='U_F')
        U[tt,...] = tmp[:,ind]
        print('\t%d of %d...'%(tt,sun.Nt))

    for mask,poly in zip(masks,polynames):
        tmp_dA  = np.sum( np.sum(U*mask,axis=-1), axis=-1)
        data[poly]['Q_r']=tmp_dA


    return data


## Testing
#################################
#ncfile = 'data/Heatflux_AVG.0'
#cellindex=range(0,9)
##cellindex=[4]
#trange = [0,-1]
####
#
#budget = calc_avg_budget(ncfile,trange,cellindex,plot=True)



