"""
VTK Unstructured Grid plotting class
"""
from .uplot import Plot

from tvtk.api import tvtk
from mayavi import mlab
import numpy as np


class PlotVTK(Plot):
    """

    """
    def __init__(self, xp, yp, cells, **kwargs):

        Plot.__init__(self, xp, yp, cells, **kwargs)

        self.maxfaces = self.cells.shape[1]

        self.points = np.column_stack((self.xp,self.yp,0.0*self.xp))

        if self.maxfaces == 3:
            self._ug = self._init_tri_2D()
        else:
            self._ug = self._init_mixed_2D()

    ###
    # Plotting stuff
    def plotcelldata(self, data, vmin=None, vmax=None, **kwargs):
        """
        Surface plot of the scalar in the 'data' attribute with mayavi
        
        Works on the 2D and 3D data
        """

        self._ug.cell_data.scalars = data
        self._ug.cell_data.scalars.name = 'suntans_scalar'

        # Create a new scene if there isn't one
        if 'fig' not in self.__dict__:
            self.newscene()
        
        src = mlab.pipeline.add_dataset(self._ug)
        self.h=mlab.pipeline.surface(src, vmin=vmin, vmax=vmax, **kwargs)
        
        # Add a colorbar if the isn't one
        if 'cb' not in self.__dict__:
            self.colorbar() 
            
        ## Add a title if there isn't one
        #if not self.__dict__.has_key('title'):
        #    self.title=mlab.title(Spatial.genTitle(self),height=0.95,size=0.15)

        return src
        
    def contour(self,vv=[10],clim=None,**kwargs):
        """
        Filled contour plot of scalar data
        """
        
        if clim is None:
            clim = [data.min(), data.max()]
        
        # Create a new scene if there isn't one
        if 'fig' not in self.__dict__:
            self.newscene()
        
        # Need to set use point (vertex) data
        src = mlab.pipeline.cell_to_point_data(self._ug)
        
        # Add the contour_surface module to the scene
        self.h=mlab.pipeline.contour_surface(src, contours=vv, line_width=1.0, \
                vmax=clim[1],vmin=clim[0],**kwargs)
        self.h.contour.filled_contours=True # This is the trick to fill the contours
        
        # Add a colorbar if the isn't one
        if 'cb' not in self.__dict__:
            self.colorbar() 
            
        # Add a title if there isn't one
        if 'title' not in self.__dict__:
            self.title=mlab.title(Spatial.genTitle(self),height=0.95,size=0.15)

    ####
    # Initialisation functions
    def _init_tri_2D(self):
        """
        Initialise the actual 2 dimensional tvtk object
        """
        
        tri_type = tvtk.Triangle().cell_type
        
        ug = tvtk.UnstructuredGrid(points=self.points)
        ug.set_cells(tri_type, self.cells)
    
        #ug.cell_data.scalars = self.data
        #ug.cell_data.scalars.name = 'suntans_scalar'
        
        return ug
 

    def _init_mixed_2D(self):
        """
        Initialise the actual 2 dimensional tvtk object

        This is for a mixed data type
        """

        poly_type = tvtk.Polygon().cell_type
        
        ug = tvtk.UnstructuredGrid(points=self.points)
        
        # Fill all of the cells with the first points
        #    this is a bit of a hack but works...
        cells = self.cells
        for ii in range(self.Nc):
            nf = self.nfaces[ii]
            cells[ii,nf::] = cells[ii,0]
            
        #offsets, cells = self.to_metismesh()
        ##cells = np.array(self.cells[self.cells.mask==False])
        #cell_array = tvtk.CellArray()
        #cell_array.set_cells(self.Nc,cells)
        ##offsets = np.cumsum(self.nfaces)
        ##offsets = offsets - offsets[0]
        #poly_types = [poly_type for ii in range(self.Nc)]
        ###
        #ug.set_cells(np.array(poly_types), offsets, cell_array)
        ### For a single cell type
        
        ug.set_cells(poly_type, cells)
    
        #ug.cell_data.scalars = self.data
        #ug.cell_data.scalars.name = 'suntans_scalar'
        
        return ug

    
    def newscene(self,size=(800,700)):
        """
        Creates a new scene
        """
        
        self.fig=mlab.figure(bgcolor=(0.,0.,0.),size=size)
        self.fig.scene.z_plus_view()
        
    def colorbar(self):
        """
        Adds a colorbar for the object in 'h'
        """
        self.cb = mlab.colorbar(object=self.h,orientation='vertical')
    

