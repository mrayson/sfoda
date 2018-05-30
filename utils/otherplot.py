"""
Other plotting routines outside of matplotlib
"""
import matplotlib.transforms as transforms
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib import collections
from matplotlib.ticker import Formatter

from datetime import datetime, timedelta
import numpy as np

from soda.utils import othertime

class StreakPlot(object):
    """
    Class for generating a streak plot
    """

    # Width of the start and end of the tail
    widthmin = 0.1
    widthmax = 1.5
    # Plotting parameters
    colors = '0.9'
    alpha = 0.9

    xlim=None
    ylim=None

    def __init__(self,xp,yp,**kwargs):
        """
        StreakPlot:

        xp, yp - spatial location matrices with dimension Nt x Nxy

        Streaks are along the first dimension Nt
        """

        self.__dict__.update(**kwargs)

        self.xp = xp
        self.yp = yp
        self.nx,self.nt = np.shape(self.xp)

        self.create_xy()


    def create_xy(self):
        # Create the start and end points
        self.points = np.concatenate([self.xp[...,np.newaxis],\
            self.yp[...,np.newaxis]],axis=-1)

        # Create the linewidths
        lwidths = np.linspace(self.widthmin,self.widthmax,self.nt-1)
        lwidths = lwidths**2 # this thins out the tail a bit
        self.lwidths = lwidths/lwidths[-1] * self.widthmax

        # Create a list of segment arrays
        self.segments=[ np.concatenate([self.points[ii,:-1,np.newaxis,:],\
            self.points[ii,1:,np.newaxis,:]],axis=1)\
            for ii in range(self.nx) ]

        # Create a list of line collections
        self.lcs = [ LineCollection(segment, linewidths=self.lwidths,\
            colors=self.colors,alpha=self.alpha)\
            for segment in self.segments ]

    def update(self,xp,yp):
        """
        Update the line collections with new x,y points
        """
        self.xp=xp
        self.yp=yp
        # Create the start and end points
        self.points = np.concatenate([xp[...,np.newaxis],yp[...,np.newaxis]],axis=-1)

        # Create a list of segment arrays
        self.segments=[ np.concatenate([self.points[ii,:-1,np.newaxis,:],\
            self.points[ii,1:,np.newaxis,:]],axis=1)\
            for ii in range(self.nx) ]

        for lc, seg in zip(self.lcs,self.segments):
            lc.set_segments(seg)

    def plot(self, ax):
        """Inserts each line collection into current plot"""

        if self.xlim == None:
            self.xlim = [self.xp.min(),self.xp.max()]
        if self.ylim == None:
            self.ylim = [self.yp.min(),self.yp.max()]

        ax.set_xlim(self.xlim)
        ax.set_ylim(self.ylim)
        ax.set_aspect('equal')

        list(map(ax.add_collection,self.lcs))
        #for lc in self.lcs:
        #    ax.add_collection(lc)

        return ax

def streakplot(xp,yp,ax=None,**kwargs):
    """
    Functional call to StreakPlot class
    """
    if ax==None:
        ax=plt.gca()

    S = StreakPlot(xp,yp,**kwargs)
    S.plot(ax)
    return S

def stackplot(t,y,scale=None,gap=0.2,ax=None,fig=None,units='',labels=None,**kwargs):
    """
    Vertically stacked time series plot.
    
    Puts all of the time-series into one axes by working out a suitable spacing. 
    
    Inputs:
        y - 2d array [nt,ny] where ny is the number of time series
        t - datetime vector
        
    Returns: 
        fig, ax : figure and axes handles
        ll : plot handles to each line plot [list]
    """
    # Determine the scale factors and the heights of all of the axes
    ny = y.shape[0]

    # Make sure that the time is datetime
    if isinstance(t[0], np.datetime64):
        t = othertime.datetime64todatetime(t)
        
    if scale==None:
        scale = np.abs(y).max()
    
    if not labels == None:
        assert len(labels)==ny, ' number of labels (%d) must equal number of layers (%d)'%(len(labels),ny)
        
    # Height of each axes in normalized coordinates
    yheight = 1.0 / (ny + (ny+1.0)*gap)
    
    # Create a new figure
    if fig==None:
        fig=plt.figure()
    else:
        fig = plt.gcf()
    
    if ax == None:
        ax = fig.add_subplot(111,frame_on=False,ylim=[0,1.0],yticks=[])
        
    # Now add each line to the figure
    ll = [] # List of line objects
    
    def fakeaxes(yval,dy):
        cc=[0.5,0.5,0.5]
        ax.add_line(Line2D([0,1],[yval,yval],linewidth=0.5,color=cc,transform=ax.transAxes,linestyle='--'))
        yp = yval + dy/2.
        ym = yval - dy/2.
        ax.add_line(Line2D([0,0],[yp,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([1,1],[yp,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        #Little caps
        ax.add_line(Line2D([0,0.01],[yp,yp],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0,0.01],[ym,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0.99,1],[yp,yp],linewidth=0.5,color=cc,transform=ax.transAxes))
        ax.add_line(Line2D([0.99,1],[ym,ym],linewidth=0.5,color=cc,transform=ax.transAxes))
        
    for N in range(1,ny+1):
        yoffset = N*(gap*yheight) + 0.5*yheight + (N-1)*yheight     
        # scaling factor
        #vscale = yheight / (scale+yoffset)
        vscale = yheight / (2*scale)
        l = ax.plot(t,vscale*y[N-1,:]+yoffset,**kwargs)
        ll.append(l)
        #Adds an axes
        fakeaxes(yoffset,yheight)
        
        if not labels==None:
            plt.text(0.2,yoffset+0.5*yheight-0.02,labels[N-1],transform=ax.transAxes,fontstyle='italic')
              
    # Add a few extra features    
    ax.add_line(Line2D([0,1],[0.01,0.01],linewidth=0.5,color='k',transform=ax.transAxes))
    ax.add_line(Line2D([0,1],[1,1],linewidth=0.5,color='k',transform=ax.transAxes))
    plt.xticks(rotation=17)
    plt.ylabel('Scale = $\pm$%2.1f  [%s]'%(scale,units))
    
    return fig,ax,ll

def polar_pdf( u, v,\
        speedmax=1.0, ndirbins=90, nspeedbins=25, cmap='RdYlBu_r',\
        ax=None):
        """
        Polar plot of speed-direction joint frequency distribution
        """
        # Convert cartesian polar coordinates
        u_hat_mod = u + 1j*v
        speed_mod = np.abs(u_hat_mod)
        theta_mod = np.angle(u_hat_mod)
 
        # Create a 2d histogram of speed direction
        abins = np.linspace(-np.pi, np.pi, ndirbins)      # 0 to 360 in steps of 360/N.
        sbins = np.linspace(0.0, speedmax, nspeedbins) 

        Hmod, xedges, yedges = np.histogram2d(theta_mod, speed_mod,\
                bins=(abins,sbins), normed=True)

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

        #fig = plt.figure( figsize = (14,6))

        if ax is None:
            ax = plt.subplot(111, projection='polar')

        ax.set_xticklabels(xlabels)
        C = ax.contourf(theta, r, Hmod.T, clevs, cmap = cmap)

        return C, ax

    
def SpeedDirPlot(t,u,v,convention='current',units='m s^{-1}',color1='b',color2='r'):
    """
    Plots speed and direction on the same axes
    
    Inputs:
        t - time vector
        u,v - velocity cartesian components
        
    Returns:
        ax - list of axes  handles
        h - list of plot handles
    
    convention = 'current' or 'wind'
    
    See this example:
        http://matplotlib.org/examples/api/two_scales.html
    """
    import airsea
    
    Dir, Spd = airsea.convertUV2SpeedDirn(u,v,convention=convention)
    
    
    ax = list(range(2))
    h = list(range(2))
    fig = plt.gcf()
    ax[0] = fig.gca()
    
    
    # Left axes
    h[0] = ax[0].fill_between(t, Spd, color=color1,alpha=0.7)
    # Make the y-axis label and tick labels match the line color.
    ax[0].set_ylabel('Speed [$%s$]'%units, color=color1)
    for tl in ax[0].get_yticklabels():
        tl.set_color(color1)

    #Right axes
    ax[1] = ax[0].twinx() # This sets up the second axes
    ax[1].plot(t, Dir, '.',color=color2)
    ax[1].set_ylabel("Dir'n [$\circ$]", color=color2)
    ax[1].set_ylim([0,360])
    ax[1].set_yticks([0,90,180,270])
    ax[1].set_yticklabels(['N','E','S','W'])
    for tl in ax[1].get_yticklabels():
        tl.set_color(color2)
        
    plt.setp( ax[0].xaxis.get_majorticklabels(), rotation=17 )
        
    return ax, h

def ProfilePlot(t,y,z,scale=86400,\
        axis=0,color=[0.5,0.5,0.5],xlim=None,units='m/s',scalebar=1.0):
    """
    Plot a series of vertical profiles as a time series
    
    scale - Sets 1 unit = scale (seconds)
    
    See this page on formatting:
        http://matplotlib.org/examples/pylab_examples/date_index_formatter.html
    """
    
    class MyFormatter(Formatter):
        def __init__(self, dates, fmt='%b %d %Y'):
            self.fmt = fmt
            self.dates = dates

        def __call__(self, x, pos=0):
            'Return the label for time x s'
            return datetime.strftime(datetime(1990,1,1)+timedelta(seconds=x),self.fmt)

    tsec = othertime.SecondsSince(t)
    formatter = MyFormatter(tsec)
    
    y = np.swapaxes(y,0,axis)
    
    lines=[]
    line2 =[]
    for ii, tt in enumerate(tsec):
        #xplot = set_scale(y[:,ii],tt)
        xplot = tt + y[:,ii]*scale
        lines.append(np.array((xplot,z)).T)
        line2.append(np.array([[tt,tt],[z[0],z[-1]]]).T)
        
    
    LC1 = collections.LineCollection(lines,colors=color,linewidths=1.5)
    LC2 = collections.LineCollection(line2,colors='k',linestyles='dashed') # Zero axis
    
    ax=plt.gca()
    ax.add_collection(LC1)
    ax.add_collection(LC2)
    ax.set_ylim((z.min(),z.max()))
    ax.xaxis.set_major_formatter(formatter)
    if xlim==None:
        xlim=(tsec[0]-scale/2,tsec[-1]+scale/2)
    else:
        xlim=othertime.SecondsSince(xlim)
    ax.set_xlim(xlim)
    plt.xticks(rotation=17)       

    ###
    # Add a scale bar    
    ###
    
    # Compute the scale bar size in dimensionless units
    if not scalebar==None:
        xscale = scalebar*scale/(xlim[-1]-xlim[0])
        x0 = 0.1
        y0 = 0.8
        dy = 0.02
        ax.add_line(Line2D([x0,x0+xscale],[y0,y0],linewidth=0.5,color='k',transform=ax.transAxes))
        #Little caps
        ax.add_line(Line2D([x0,x0],[y0-dy,y0+dy],linewidth=0.5,color='k',transform=ax.transAxes))
        ax.add_line(Line2D([x0+xscale,x0+xscale],[y0-dy,y0+dy],linewidth=0.5,color='k',transform=ax.transAxes))
        plt.text(x0,y0+0.05,'Scale %3.1f %s'%(scalebar,units),\
            transform=ax.transAxes)

    
    return ax
    
def monthlyhist(t,y,ylim=0.1,xlabel='',ylabel='',title='',**kwargs):
    """
    Plots 12 histograms on a 6x2 matrix of variable, y, grouped by calendar month
    
    Inputs:
        y - vector of data
        t - vector of datetime objects
        kwargs - keyword arguments for numpy.hist
            
    """
    month = othertime.getMonth(t)
    fig=plt.gcf()
    
    for m in range(1,13):
        
        # Find the values
        ind = np.argwhere(month==m)
        data=y[ind]
        
        ax=plt.subplot(6,2,m)
        if len(data)>0:
            plt.hist(data,**kwargs)
        
        mon=datetime.strftime(datetime(1900,m,1),'%B')
        plt.title(mon)
        plt.ylim([0,ylim]) 
        
    
        if m not in (11,12):         
            ax.set_xticklabels([])
        else:
            plt.xlabel(xlabel)
            
            
        if m not in (1,3,5,7,9,11):
            ax.set_yticklabels([])
        else:
            plt.ylabel(ylabel)
            
        
        #Calc some stats
        textstr = 'Mean: %6.1f\nStd. Dev.: %6.1f\n'%(np.mean(data),np.std(data))
        plt.text(0.5,0.5,textstr,transform=ax.transAxes)
        
        # plot a title
        plt.figtext(0.5,0.95,title,fontsize=14,horizontalalignment='center')
        
    return fig
    
def window_index(serieslength,windowsize,overlap):
    """
    Determines the indices for start and end points of a time series window
    
    Inputs:
        serieslength - length of the vector [int]
        windowsize - length of the window [int]
        overlap - number of overlap points [int]
        
    Returns: pt1,pt2 the start and end indices of each window
    """

    p1=0
    p2=p1 + windowsize
    pt1=[p1]
    pt2=[p2]
    while p2 < serieslength:
        p1 = p2 - overlap
        p2 = min((p1 + windowsize, serieslength))
        pt1.append(p1)
        pt2.append(p2)
        
    return pt1, pt2
    


def axcolorbar(cbobj,pos=[0.7, 0.8, 0.2, 0.04],ax=None,fig=None,orientation='horizontal',**kwargs):
	"""
	Inserts a colorbar with a position relative to an axes and not a figure
	
	Inputs:
		cbobj - plot object for colorbar
		pos - position vector [x0, y0, width, height] in dimensionless coordinates
		ax - axes to insert colobar
		figure - figure 
		**kwargs - arguments for plt.colorbar
	
	Returns a colorbar object
	
	Derived from this post:
		http://stackoverflow.com/questions/22413211/cant-fix-position-of-colorbar-in-image-with-multiple-subplots
	"""
	if fig is None:
		fig=plt.gcf()
	if ax is None:
		ax=plt.gca()
		
	fig.tight_layout()  # You call fig.tight_layout BEFORE creating the colorbar

	# You input the POSITION AND DIMENSIONS RELATIVE TO THE AXES
	x0, y0, width, height = pos

	# and transform them after to get the ABSOLUTE POSITION AND DIMENSIONS
	Bbox = transforms.Bbox.from_bounds(x0, y0, width, height)
	trans = ax.transAxes + fig.transFigure.inverted()
	l, b, w, h = transforms.TransformedBbox(Bbox, trans).bounds

	# Now just create the axes and the colorbar
	cbaxes = fig.add_axes([l, b, w, h])
	cbar = plt.colorbar(cbobj, cax=cbaxes,orientation=orientation, **kwargs)
	cbar.ax.tick_params(labelsize=9)

	return cbar
