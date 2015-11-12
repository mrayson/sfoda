"""
Unstructured grid generation utilities
"""
import numpy as np
from soda.dataio.ugrid.hybridgrid import HybridGrid
from soda.utils.inpolygon import inpolygon


def cartesian_ugrid_gen(xlims, ylims, dx, suntanspath=None, maskpoly=None):
    """
    Creates a cartesian grid in a box 

    Inputs:
        xlims - [xmin, xmax] 
        ylims - [ymin, ymax] 
        dx - grid spacing
    """

    #xgrd = np.arange(xlims[0]-dx/2,xlims[1]+1.5*dx,dx)
    #ygrd = np.arange(ylims[0]-dx/2,ylims[1]+1.5*dx,dx)
    xgrd = np.arange(xlims[0],xlims[1]+1.0*dx,dx)
    ygrd = np.arange(ylims[0],ylims[1]+1.0*dx,dx)

    nx = xgrd.shape[0]
    ny = ygrd.shape[0]

    # Create a mask polygon
    X,Y = np.meshgrid(xgrd,ygrd)

    return curv_ugrid_gen(X,Y,suntanspath=suntanspath,maskpoly=maskpoly)


def curv_ugrid_gen(X,Y,suntanspath=None,maskpoly=None):
    """
    Creates a curvilinear mesh from grid corners points stored
    in arrays X and Y.
    """

    ny,nx = X.shape

    XY = np.vstack((X.ravel(),Y.ravel())).T
    if not maskpoly is None:
        mask = inpolygon(XY,maskpoly)
        mask = mask.reshape((ny,nx))
    else:
        mask = np.ones((ny,nx),dtype=np.bool) # all false

    cells=[]
    xp = []
    yp = []

    def pntindx(j,i,nrows):
        return j*nrows + i

    for jj in range(ny):
        for ii in range(nx):
            #if mask[jj,ii]:
            #xp.append(xgrd[ii])
            #yp.append(ygrd[jj])
            xp.append(X[jj,ii])
            yp.append(Y[jj,ii])

    for jj in range(ny-1):
        for ii in range(nx-1):
            if mask[jj,ii] and mask[jj+1,ii] and mask[jj+1,ii+1] and mask[jj,ii+1]:
                cells.append([pntindx(jj,ii,nx), pntindx(jj+1,ii,nx),\
                    pntindx(jj+1,ii+1,nx),pntindx(jj,ii+1,nx)])
            
    Nc = len(cells)
    nfaces = 4*np.ones((Nc,),np.int)

    cells = np.array(cells,dtype=np.int32)
    # Convert to a suntans grid
    grd = HybridGrid(xp,yp,cells,nfaces=nfaces)

    if not suntanspath is None:
        grd.write2suntans(suntanspath)

    return grd
