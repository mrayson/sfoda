#!/usr/bin/python
"""
SUNTANS NetCDF plotting GUI
"""
import os
#import wx

# The recommended way to use wx with mpl is with the WXAgg
# backend. 
#
import matplotlib

#matplotlib.use('WXAgg')
# Set some default parameters
matplotlib.rcParams['text.color']='white'
matplotlib.rcParams['savefig.facecolor']='black'
matplotlib.rcParams['savefig.edgecolor']='black'
matplotlib.rcParams['figure.facecolor']='black'
matplotlib.rcParams['figure.edgecolor']='black'
matplotlib.rcParams['axes.facecolor']='black'
matplotlib.rcParams['axes.edgecolor']='white'
matplotlib.rcParams['axes.labelcolor']='white'
matplotlib.rcParams['xtick.color']='white'
matplotlib.rcParams['ytick.color']='white'
matplotlib.rcParams['font.family']='serif'
#matplotlib.rcParams['font.sans-serif']='Arial'

from matplotlib.figure import Figure
from matplotlib.collections import PolyCollection, LineCollection
import matplotlib.animation as animation
#from matplotlib.backends.backend_wxagg import \
#    FigureCanvasWxAgg as FigCanvas, \
#    NavigationToolbar2WxAgg as NavigationToolbar

from PyQt5 import QtGui
from PyQt5.QtWidgets import QAction,\
        QFileDialog, QWidget,\
        QSlider, QMainWindow,\
        QApplication, QVBoxLayout,\
        QHBoxLayout, QComboBox, QCheckBox, QLineEdit

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.uic import loadUiType

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import matplotlib

import sys


from soda.dataio.suntans.sunpy import Spatial, Grid
from soda.dataio.suntans.sunxray import Sunxray, Sundask
from datetime import datetime
import numpy as np

# Dask stuff
#from dask.distributed import Client, LocalCluster, worker, progress, wait
#
## Try and connect to a scheduler previously launched
#print('Connecting to dask scheduler...')
#client = Client(scheduler_file="~/scheduler.json", asynchronous=False)

# Use cmocean colormaps (these are awesome!!)
#try:
#    from cmocean import cm
#    USECMOCEAN=True
#except:
#    print 'No cmocean'
#    USECMOCEAN=False
USECMOCEAN=False

import pdb

#class SunPlotPy(QMainWindow, Spatial, Grid ):
class SunPlotPyX(Sundask, QMainWindow):
    """ 
    The main frame of the application
    """
    title = 'sunplot(py)'

    # Plotting options
    autoclim=True
    showedges=False
    bgcolor='k'
    textcolor='w'
    cmap='RdBu'
    particlesize = 1.8
    particlecolor = 'm'

    # other flags
    collectiontype='cells'
    oldcollectiontype='cells'

    # 
    tstep=0 
    tstepold = 0
    klayer = [0]
    variable = 'dv'
    depthlevs = [0., 10., 100., 200., 300., 400., 500.,\
        1000.,2000.,3000.,4000.,5000]

    _FillValue=999999
    
    def __init__(self, parent=None):
        #wx.Frame.__init__(self, None, -1, self.title)
        QMainWindow.__init__(self, parent)
        #super(SunPlotPy, self).__init__(parent)
        
        self.create_menu()
        #self.create_status_bar()
        self.create_main_panel()
        
        #self.draw_figure()


    def create_menu(self):
        self.file_menu = self.menuBar().addMenu("&File")
        
        ###
        # File Menu
        ###
        # Load a hydro output file
        load_file_action = self.create_action("&Open file",\
                shortcut="ctrl-o", slot=self.on_open_file,
                tip="open netcdf file")

        load_pfile_action = self.create_action("&Open parallel files",\
                slot=self.on_open_p_file,
                tip="open netcdf file")


        load_grid_action = self.create_action("&Open grid",\
                shortcut="ctrl-g", slot=self.on_load_grid,
                tip="open suntans grid")

        save_anim_action = self.create_action("&Save animation",\
                shortcut="ctrl-a", slot=self.on_save_anim,
                tip="animate current scene")

        quit_action = self.create_action("&Exit",\
                shortcut="ctrl-x", slot=self.close,
                tip="Close the application")

        self.add_actions(self.file_menu,
                (load_file_action, load_pfile_action, load_grid_action,\
                save_anim_action, None, quit_action))


        #    self.Bind(wx.EVT_MENU, self.on_open_file, m_expt)
        ###
        # Tools menu
        ###
        self.tools_menu = self.menuBar().addMenu("&Tools")
        
        ###
        # File Menu
        ###
        # Load a hydro output file
        load_stat_action = self.create_action("&Plot grid size statistics",\
                slot=self.on_plot_gridstat,
                tip="grid stats")

        self.add_actions(self.tools_menu,
                (load_stat_action,))





    #    # Load a grid file
    #    m_grid = menu_file.Append(-1, "&Load grid\tCtrl-G", "Load SUNTANS grid from folder")
    #    self.Bind(wx.EVT_MENU, self.on_load_grid, m_grid)

    #    # Load a particle file
    #    m_part = menu_file.Append(-1, "&Load PTM file\tCtrl-Shift-P", "Load a PTM file")
    #    self.Bind(wx.EVT_MENU, self.on_load_ptm, m_part)

    #    # Save current scene as an animation
    #    m_anim = menu_file.Append(-1,"&Save animation of current scene\tCtrl-S","Save animation")
    #    self.Bind(wx.EVT_MENU, self.on_save_anim, m_anim)

    #    # Save the current figure
    #    m_prin = menu_file.Append(-1,"&Print current scene\tCtrl-P","Save figure")
    #    self.Bind(wx.EVT_MENU, self.on_save_fig, m_prin)



    #    menu_file.AppendSeparator()
    #    # Exit
    #    m_exit = menu_file.Append(-1, "E&xit\tCtrl-X", "Exit")
    #    self.Bind(wx.EVT_MENU, self.on_exit, m_exit)

    #    ###
    #    # Tools menu
    #    ###
    #    menu_tools = wx.Menu()
    #    m_gridstat = menu_tools.Append(-1, "&Plot grid size statistics", "SUNTANS grid size")
    #    self.Bind(wx.EVT_MENU, self.on_plot_gridstat, m_gridstat)

    #    m_countcells = menu_tools.Append(-1, "&Count # grid cells", "Grid cell count")
    #    self.Bind(wx.EVT_MENU, self.on_count_cells, m_countcells)

    #    m_overlaybathy = menu_tools.Append(-1, "&Overlay depth contours", "Depth overlay")
    #    self.Bind(wx.EVT_MENU, self.on_overlay_bathy, m_overlaybathy)

    #    
    #    ###
    #    # Help Menu
    #    ###
    #    menu_help = wx.Menu()
    #    m_about = menu_help.Append(-1, "&About\tF1", "About the demo")
    #    self.Bind(wx.EVT_MENU, self.on_about, m_about)
    #    
    #    
    #    # Add all of the menu bars
    #    self.menubar.Append(menu_file, "&File")
    #    self.menubar.Append(menu_tools, "&Tools")
    #    self.menubar.Append(menu_help, "&Help")
    #    self.SetMenuBar(self.menubar)

    def add_actions(self, target, actions):
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)


    def create_action(  self, text, slot=None, shortcut=None, 
                    icon=None, tip=None, checkable=False, 
                    ):
                    #signal="triggered()"):
        action = QAction(text, self)
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % icon))
        if shortcut is not None:
            action.setShortcut(shortcut)

        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            #self.connect(action, SIGNAL(signal), slot)
            # Qt5 style
            action.triggered.connect(slot)
        if checkable:
            action.setCheckable(True)
        return action



    def create_main_panel(self):
        """ Creates the main panel with all the controls on it:
             * mpl canvas 
             * mpl navigation toolbar
             * Control panel for interaction
        """
        self.panel = QWidget()
        
        # Create the mpl Figure and FigCanvas objects. 
        # 5x4 inches, 100 dots-per-inch
        #
        self.dpi = 100
        #self.fig = Figure((7.0, 6.0), dpi=self.dpi,facecolor=self.bgcolor)
        self.fig = Figure((7.0, 6.0), dpi=self.dpi)
        #self.canvas = FigCanvas(self.panel, -1, self.fig)
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setParent(self.panel)
        
        # Since we have only one plot, we can use add_axes 
        # instead of add_subplot, but then the subplot
        # configuration tool in the navigation toolbar wouldn't
        # work.
        #
        self.axes = self.fig.add_subplot(111)
        #SetAxColor(self.axes,self.textcolor,self.bgcolor)
        
        # Bind the 'pick' event for clicking on one of the bars
        #
        #self.canvas.mpl_connect('pick_event', self.on_pick)
        
        #########
        ## Create widgets
        #########
        self.variable_list = QComboBox()
        self.variable_list.addItem("Select a variable...")
        self.variable_list.activated[str].connect(self.on_select_variable)

        self.time_list = QComboBox()
        self.time_list.addItem("Select a time...")
        self.time_list.activated[int].connect(self.on_select_time)

        self.depthlayer_list = QComboBox()
        self.depthlayer_list.addItem("Select a vertical layer...")
        self.depthlayer_list.activated[int].connect(self.on_select_depth)

        self.show_edge_check = QCheckBox('Show Edges', self)
        #self.show_edge_check.toggle()
        self.show_edge_check.stateChanged.connect(self.on_show_edges)

        cmaps = list(matplotlib.cm.datad.keys())
        cmaps.sort()
        self.colormap_list = QComboBox()
        self.colormap_list.clear()
        self.colormap_list.addItems(cmaps)
        self.colormap_list.activated[str].connect(self.on_select_cmap)


        self.clim_check = QCheckBox('Manual color limits', self)
        #self.show_edge_check.toggle()
        self.clim_check.stateChanged.connect(self.on_clim_check)

        #self.clim_check = wx.CheckBox(self.panel, -1, 
        #    "Manual color limits ",
        #    style=wx.ALIGN_RIGHT)
        #self.clim_check.Bind(wx.EVT_CHECKBOX, self.on_clim_check)

        self.climlow = QLineEdit()
        self.climlow.textChanged[str].connect(self.on_climlow)

        self.climhigh = QLineEdit()
        self.climhigh.textChanged[str].connect(self.on_climhigh)


        ## Labels
        #self.variable_label = wx.StaticText(self.panel, -1,"Variable:",size=(200,-1))
        #self.time_label = wx.StaticText(self.panel, -1,"Time step:",size=(200,-1))
        #self.depth_label = wx.StaticText(self.panel, -1,"Vertical level:",size=(200,-1))


        # Create the navigation toolbar, tied to the canvas
        #
        self.toolbar = NavigationToolbar(self.canvas, self)
        #self.toolbar.toolitems[8][3]='my_save_fig'

        #def my_save_fig(self,*args):
        #     print 'saving figure'
        ##    return "break"

        
        #########
        # Layout with box sizers
        #########
        hbox = QHBoxLayout()
        for w in [self.variable_list, self.time_list, self.depthlayer_list]:
            hbox.addWidget(w)
            hbox.setAlignment(w, Qt.AlignVCenter)

        hbox1 = QHBoxLayout()
        for w in [self.show_edge_check, self.colormap_list, self.clim_check,
                self.climlow, self.climhigh]:
            hbox1.addWidget(w)
            hbox1.setAlignment(w, Qt.AlignVCenter)


        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.canvas)
        self.vbox.addWidget(self.toolbar)
        self.vbox.addLayout(hbox)
        self.vbox.addLayout(hbox1)

        self.panel.setLayout(self.vbox)
        self.setCentralWidget(self.panel)

   #
    ###########
    ## Event functions
    ###########

    def gen_title(self,tt=None):
        
        if tt  is None:
            if type(self.tstep)==int:
                tt = self.tstep
            else:
                tt = self.tstep[0]
            
        if self.klayer[0] in [-1,'seabed']:
            zlayer = 'seabed'
        elif self.klayer[0] =='surface':
            zlayer = 'surface'
        elif self.klayer[0]==-99:
            zlayer = 'depth-averaged'
        elif self.klayer[0]>=0:
            zlayer = '%3.1f [m]'%self._ds.z_r[self.klayer[0]]

        if 'time' in self.__dict__:
            #tstr = datetime.strftime(self._ds.time[tt],\
            #    '%Y-%m-%d %H:%M:%S')
            str = self._ds.time.values[tt].astype(str)[0:-10]
            tstr = 'Time: %s'%tstr
        else:
            tstr = ''
        var = self._ds[self.variable]
        titlestr='%s [%s]\n z: %s, %s'%(var.long_name, var.units, zlayer, tstr)
                
        return titlestr



    def create_figure(self):
        """ 
        Creates the figure
        """
        # Find the colorbar limits if unspecified
        if self.autoclim:
            self.clim = [self.data.min(),self.data.max()]
            self.climlow.setText('%3.1f'%self.clim[0])
            self.climhigh.setText('%3.1f'%self.clim[1])

         
        if 'collection' in self.__dict__:
            #self.collection.remove()
            self.axes.collections.remove(self.collection)
        else:
            # First call - set the axes limits
            self.axes.set_aspect('equal')
            self.axes.set_xlim(self.xlims)
            self.axes.set_ylim(self.ylims)
 

        if self.collectiontype=='cells':
            self.collection = PolyCollection(self.xy(),cmap=self.cmap)
            self.collection.set_array(np.array(self.data[:]))
            if not self.showedges:
                self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
        elif self.collectiontype=='edges':
            xylines = [self.xp[self.edges],self.yp[self.edges]]
            linesc = [list(zip(xylines[0][ii,:],xylines[1][ii,:])) for ii in range(self.Ne)]
            self.collection = LineCollection(linesc,array=np.array(self.data[:]),cmap=self.cmap)

        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])

        self.axes.add_collection(self.collection)    
        self.title=self.axes.set_title(self.gen_title(),color=self.textcolor)
        self.axes.set_xlabel('Easting [m]')
        self.axes.set_ylabel('Northing [m]')

        # create a colorbar

        if 'cbar' not in self.__dict__:
            self.cbar = self.fig.colorbar(self.collection)
            #SetAxColor(self.cbar.ax.axes,self.textcolor,self.bgcolor)
        else:
            #pass
            print('Updating colorbar...')
            #self.cbar.check_update(self.collection)
            self.cbar.on_mappable_changed(self.collection)

        self.canvas.draw()
   
    def update_figure(self):
        if self.autoclim:
            self.clim = [self.data.min(),self.data.max()]
            self.climlow.setText('%3.1f'%self.clim[0])
            self.climhigh.setText('%3.1f'%self.clim[1])
        else:
            self.clim = [float(self.climlow.text()),\
                float(self.climhigh.text())]
 
        # check whether it is cell or edge type
        if self.has_dim(self.variable,'Ne'):
            self.collectiontype='edges'
        elif self.has_dim(self.variable,'Nc'):
            self.collectiontype='cells'

        # Create a new figure if the variable has gone from cell to edge of vice
        # versa
        if not self.collectiontype==self.oldcollectiontype:
            self.create_figure()
            self.oldcollectiontype=self.collectiontype

        self.collection.set_array(np.array(self.data[:]))
        self.collection.set_clim(vmin=self.clim[0],vmax=self.clim[1])

        # Cells only
        if self.collectiontype=='cells':
            if not self.showedges:
                self.collection.set_edgecolors(self.collection.to_rgba(np.array((self.data[:])))) 
                #self.collection.set_edgecolors(None)
            else:
                self.collection.set_edgecolors('k')
                self.collection.set_linewidths(0.2)

        # Update the title
        self.title=self.axes.set_title(self.gen_title(),color=self.textcolor)

        #Update the colorbar
        self.cbar.update_normal(self.collection)

        # redraw the figure
        self.canvas.draw()
    
    #def on_pick(self, event):
    #    # The event received here is of the type
    #    # matplotlib.backend_bases.PickEvent
    #    #
    #    # It carries lots of information, of which we're using
    #    # only a small amount here.
    #    # 
    #    box_points = event.artist.get_bbox().get_points()
    #    msg = "You've clicked on a bar with coords:\n %s" % box_points
    #    
    #    dlg = wx.MessageDialog(
    #        self, 
    #        msg, 
    #        "Click!",
    #        wx.OK | wx.ICON_INFORMATION)

    #    dlg.ShowModal() 
    #    dlg.Destroy()        
    #
    def on_select_variable(self, event):
        #vname = event.GetString()
        vname = event
        #self.flash_status_message("Selecting variable: %s"%vname)
        # update the spatial object and load the data
        self.variable = vname
        #self.loadData(variable=self.variable)
        self.data = self.load_data(self.variable, tstep=self.tstep, klayer=self.klayer[0])
        #self.data = data.values.ravel()

        # Check if the variable has a depth coordinate
        depthstr = ['']
        # If so populate the vertical layer box
        if self.has_dim(self.variable, 'Nk'):
            depthstr = ['%3.1f'%self._ds.z_r[k] for k in range(self.Nkmax)]
            depthstr += ['surface','seabed']

        elif self.has_dim(self.variable,'Nkw'):
            depthstr = ['%3.1f'%self._ds.z_w[k] for k in range(self.Nkmax+1)]

        self.depthlayer_list.clear()
        self.depthlayer_list.addItems(depthstr)

        # Update the plot
        self.update_figure()



    def on_select_time(self, event):

        self.tstep = event#
        # Update the object time index and reload the data
        #if self.plot_type=='hydro':
        if not self.tstep==self.tstepold:
            self.tstepold = self.tstep*1
            self.data = self.load_data(self.variable, tstep=self.tstep, klayer=self.klayer[0])
            #self.data = data.values.ravel()
            #self.flash_status_message("Selecting variable: %s..."%event.GetString())

            # Update the plot
            self.update_figure()
        #elif self.plot_type=='particles':
        #    self.PTM.plot(self.tindex,ax=self.axes,\
        #        xlims=self.axes.get_xlim(),ylims=self.axes.get_ylim())
        #
        #    self.canvas.draw()


    def on_select_depth(self, event):
        kindex = event
        if not self.klayer[0]==kindex:
            # Check if its the seabed or surface value
            if kindex>=self.Nkmax:
                kindex=event.GetString()
            self.klayer = [kindex]
            self.data = self.load_data(self.variable, tstep=self.tstep, klayer=self.klayer[0])
            #self.data = data.values.ravel()
            #self.flash_status_message("Selecting depth: %s..."%event.GetString())

            # Update the plot
            self.update_figure()

    def on_open_p_file(self, event):
        file_choices = "SUNTANS NetCDF (*.nc.*);;All Files (*.*)"

        dlg = QFileDialog.getOpenFileNames(self, "Open SUNTANS parallel file...",
                "", file_choices)

        path = dlg[0]

        if len(path) == 0:
            return

        fileparts = path[0].split('.')
        filecard = '%s.nc.*'%fileparts[-3]

        print(fileparts)
        print(filecard)
        self.on_open_file([filecard], newfile=False)


    def on_open_file(self, event, newfile=True):
        if newfile:
            file_choices = "SUNTANS NetCDF (*.nc);;All Files (*.*)"

            dlg = QFileDialog.getOpenFileNames(self, "Open SUNTANS file...",
                    "", file_choices)

            path = dlg[0]

            if len(path) == 0:
                return
        else:
            path = event
        


        self.statusBar().showMessage("Opening SUNTANS file: %s" % path)

        Sundask.__init__(self, path[0], )

        print('hellooo')
        self.Nkmax = self.load_data('Nk').max()-1

        self.statusBar().clearMessage()

        # Populate the drop down menus
        #vnames = self._ds.variables.keys()
        vnames = self.list_coord_vars()
        self.variable_list.clear()
        self.variable_list.addItems(vnames)
        # Update the time drop down list

        #if self.__dict__.has_key('time'):
        if 'time' in list(self._ds.variables.keys()):
            #self.timestr = [datetime.strftime(tt,'%d-%b-%Y %H:%M:%S') for tt in self._ds.time.values]
            self.timestr = [tt.astype(str)[:-10] for tt in self._ds.time.values]
        else:
            # Assume that it is a harmonic-type file
            self.timestr = self.nc.Constituent_Names.split()

        self.time_list.clear()
        self.time_list.addItems(self.timestr)

        # Draw the depth
        if self.variable in vnames:
            self.data = self.load_data(self.variable)
            self.create_figure()

    def on_load_grid(self, event):
        
        dir_ = QFileDialog.getExistingDirectory(None, 'Select a SUNTANS grid folder:',\
                '~/', QFileDialog.ShowDirsOnly)
        
        if dir_ is not None:
            path = dir_

            # Initialise the class
            #self.flash_status_message("Opening SUNTANS grid from folder: %s" % path)
            Grid.__init__(self,path)

            # Plot the Grid
            if 'collection' in self.__dict__:
                self.axes.collections.remove(self.collection)

            self.axes,self.collection = self.plotmesh(ax=self.axes,edgecolors='y')

            # redraw the figure
            self.canvas.draw()

    #def on_load_ptm(self, event):
    #    file_choices = "PTM NetCDF (*.nc)|*.nc|PTM Binary (*_bin.out)|*_bin.out|All Files (*.*)|*.*"
    #    
    #    dlg = wx.FileDialog(
    #        self, 
    #        message="Open PTM file...",
    #        defaultDir=os.getcwd(),
    #        defaultFile="",
    #        wildcard=file_choices,
    #        style= wx.FD_MULTIPLE)
    #    
    #    if dlg.ShowModal() == wx.ID_OK:
    #        self.plot_type = 'particles'
    #        path = dlg.GetPath()

    #        # Initialise the class
    #        if dlg.GetFilterIndex() == 0: #SUNTANS
    #            self.flash_status_message("Opening PTM netcdf file: %s" % path)
    #            self.PTM = PtmNC(path)
    #        elif dlg.GetFilterIndex() == 1: #PTM
    #            self.flash_status_message("Opening PTM binary file: %s" % path)
    #            self.PTM = PtmBin(path)

    #        self.Nt = self.PTM.nt
    #        
    #        # Update the time drop down list
    #        self.timestr = [datetime.strftime(tt,'%d-%b-%Y %H:%M:%S') for tt in self.PTM.time]
    #        self.time_list.SetItems(self.timestr)

    #        # Plot the first time step
    #        if self.__dict__.has_key('xlims'):
    #            self.PTM.plot(self.PTM.nt-1,ax=self.axes,xlims=self.xlims,\
    #            ylims=self.ylims,color=self.particlecolor,\
    #            fontcolor='w',markersize=self.particlesize)
    #        else:
    #            self.PTM.plot(self.PTM.nt-1,ax=self.axes,fontcolor='w',\
    #                color=self.particlecolor,markersize=self.particlesize)
    #        # redraw the figure
    #        self.canvas.draw()

    #    
    def on_show_edges(self,event):
        if event > 0:
            self.showedges = True
        else:
            self.showedges = False

        # Update the figure
        self.update_figure()

    def on_clim_check(self,event):
        if event > 0:
            self.autoclim=False
            self.update_figure()
        else:
            self.autoclim=True
       

    def on_climlow(self,event):
        try:
            self.clim[0] = float(event)
        except:
            return # do nothing

    #    self.clim[0] = event.GetString()
    #    #self.update_figure()

    def on_climhigh(self,event):
        try:
            self.clim[1] = float(event)
        except:
            return # do nothing
    #    print event
    #    self.clim[1] = event.GetString()
    #    #self.update_figure()

    def on_select_cmap(self,event):
        self.cmap=event
        self.collection.set_cmap(self.cmap)

        # Update the figure
        self.update_figure()

    #def on_save_fig(self,event):
    #    """
    #    Save a figure of the current scene to a file
    #    """
    #    file_choices = " (*.png)|*.png| (*.pdf)|*.pdf |(*.jpg)|*.jpg |(*.eps)|*eps "
    #    filters=['.png','.pdf','.png','.png']

    #    
    #    dlg = wx.FileDialog(
    #        self, 
    #        message="Save figure to file...",
    #        defaultDir=os.getcwd(),
    #        defaultFile="",
    #        wildcard=file_choices,
    #        style= wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

    #    if dlg.ShowModal() == wx.ID_OK:

    #        path = dlg.GetPath()
    #        ext = filters[dlg.GetFilterIndex()]
    #        if ext in path:
    #            outfile=path
    #        else:
    #            outfile = path+ext

    #        self.fig.savefig(outfile)

    #        


    def on_save_anim(self,event):
        """
        Save an animation of the current scene to a file
        """
        #file_choices = "Quicktime (*.mov)|*.mov| (*.gif)|*.gif| (*.avi)|*.avi |(*.mp4)|*.mp4 "
        #filters=['.mov','.gif','.avi','.mp4']
        filters = "Movie formats (*.mp4 *.avi *.gif);;All files (*.*)"

        
        dir_ = QFileDialog.getSaveFileName(None, 'Save animation to file:',\
                '~/', filters)
        
        #dlg = wx.FileDialog(
        #    self, 
        #    message="Output animation file...",
        #    defaultDir=os.getcwd(),
        #    defaultFile="",
        #    wildcard=file_choices,
        #    style= wx.FD_SAVE | wx.FD_OVERWRITE_PROMPT)

        if dir_ is not None:

            outfile = dir_[0]
            ext = outfile[-4::]
            # Create the animation
            #self.tstep = range(self.Nt) # Use all time steps for animation
            #self.animate(cbar=self.cbar,cmap=self.cmap,\
            #    xlims=self.axes.get_xlim(),ylims=self.axes.get_ylim())
            def initanim():
                return (self.title, self.collection)
                #if not self.plot_type=='particles':
                #    return (self.title, self.collection)
                #else:
                #    return (self.PTM.title,self.PTM.p_handle)

            def updateScalar(i):
                self.tstep=[i]
                self.data = self.load_data(self.variable, tstep=self.tstep, klayer=self.klayer[0])
                #self.data = data.values.ravel()
                #self.loadData()
                self.update_figure()
                return (self.title,self.collection)

                #if not self.plot_type=='particles':
                #    self.tstep=[i]
                #    self.loadData()
                #    self.update_figure()
                #    return (self.title,self.collection)
                #elif self.plot_type=='particles':
                #    self.PTM.plot(i,ax=self.axes,\
                #        xlims=self.axes.get_xlim(),ylims=self.axes.get_ylim())
                #    return (self.PTM.title,self.PTM.p_handle)

            self.anim = animation.FuncAnimation(self.fig, \
                updateScalar, init_func = initanim, frames=self.Nt, interval=50, blit=True)

            if ext=='.gif':
                self.anim.save(outfile,writer='imagemagick',fps=6)
            elif ext=='.mp4':
                print('Saving html5 video...')
                # Ensures html5 compatibility
                self.anim.save(outfile,fps=6,\
                    writer='ffmpeg',\
                    bitrate=3600,extra_args=['-vcodec','libx264'])
                    #writer='mencoder',
                    #bitrate=3600,extra_args=['-ovc','x264']) # mencoder options
            else:
                self.anim.save(outfile,writer='mencoder',fps=6,bitrate=3600)

            # Return the figure back to its status
            del self.anim
            #self.loadData()
            self.data = self.load_data(self.variable, tstep=self.tstep, klayer=self.klayer[0])
            #self.data = data.values.ravel()
            self.update_figure()
            print('Finished saving animation to %s'%outfile)
            print(72*'#')

            #if not self.plot_type=='particles':
            #    self.loadData()
            #    self.update_figure()

            # Bring up a dialog box
            #dlg2= wx.MessageDialog(self, 'Animation complete.', "Done", wx.OK)
            #dlg2.ShowModal()
            #dlg2.Destroy()

    #def on_exit(self, event):
    #    self.Destroy()
    #    
    #def on_about(self, event):
    #    msg = """ SUNTANS NetCDF visualization tool
    #    
    #        *Author: Matt Rayson
    #        *Institution: Stanford University
    #        *Created: October 2013
    #    """
    #    dlg = wx.MessageDialog(self, msg, "About", wx.OK)
    #    dlg.ShowModal()
    #    dlg.Destroy()

    #def on_count_cells(self,eveny):
    #    msg = "Total 3-D grid cells = %d"%(self.count_cells())
    #    dlg = wx.MessageDialog(self, msg, "No. cells", wx.OK)
    #    dlg.ShowModal()
    #    dlg.Destroy()

    #def on_overlay_bathy(self,event):
    #    # Plot depth contours
    #    print 'Plotting contours...'
    #    self.contourf(z=self.dv, clevs=self.depthlevs,\
    #        ax=self.axes,\
    #        filled=False, colors='0.5', linewidths=0.5, zorder=1e6)
    #    print 'Done'
   
    def on_plot_gridstat(self, event):
        """
        Plot the grid size histogram in a new figure
        """
        matplotlib.pyplot.figure()
        self.plothist()
        matplotlib.pyplot.show()


    #def create_status_bar(self):
    #    self.statusbar = self.CreateStatusBar()

    #def flash_status_message(self, msg, flash_len_ms=1500):
    #    self.statusbar.SetStatusText(msg)
    #    self.timeroff = wx.Timer(self)
    #    self.Bind(
    #        wx.EVT_TIMER, 
    #        self.on_flash_status_off, 
    #        self.timeroff)
    #    self.timeroff.Start(flash_len_ms, oneShot=True)
    #
    #def on_flash_status_off(self, event):
    #    self.statusbar.SetStatusText('')


def SetAxColor(ax,color,bgcolor):
    ax.set_axis_bgcolor(bgcolor)
    
    ax.yaxis.set_tick_params(color=color,labelcolor=color)
    ax.xaxis.set_tick_params(color=color,labelcolor=color)
    ax.yaxis.label.set_color(color)
    ax.xaxis.label.set_color(color)
    ax.spines['top'].set_color(color)
    ax.spines['bottom'].set_color(color)
    ax.spines['left'].set_color(color)
    ax.spines['right'].set_color(color)
    return ax
 

if __name__ == '__main__':
    app = QApplication(sys.argv)
    form = SunPlotPyX()
    form.show()
    app.exec_()

