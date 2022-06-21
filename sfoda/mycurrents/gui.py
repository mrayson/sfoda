#from oceanplotpy import *
from PyQt5 import QtGui
from PyQt5.QtWidgets import QAction, QFileDialog, QWidget, QSlider
from PyQt5.QtCore import Qt
from PyQt5.uic import loadUiType

from matplotlib.figure import Figure
from matplotlib.backends.backend_qt5agg import (
    FigureCanvasQTAgg as FigureCanvas,
    NavigationToolbar2QT as NavigationToolbar)
import matplotlib

import os
from datetime import datetime
from netCDF4 import Dataset
import numpy as np

import mycurrents
from mycurrents import oceanmooring

# Use cmocean colormaps (these are awesome!!)
try:
    from cmocean import cm
    USECMOCEAN=True
except:
    print('No cmocean')
    USECMOCEAN=False


# Load these from the designer interface
#installdir = '/home/suntans/code/mycurrents/'
_path = os.path.abspath(mycurrents.__file__)
_dirpath = os.path.dirname(_path)


Ui_MainWindow, QMainWindow = loadUiType('%s/gui/oceanplotpy.ui'%_dirpath)

# Load the colormaps
if USECMOCEAN:
    cmaps=[]
    for cmap in cm.cmapnames:
        cmaps.append(cmap)
        cmaps.append(cmap+'_r') # Add all reverse map options
else:
    # Use matplotlib standard
    cmaps = list(matplotlib.cm.datad.keys())
cmaps.sort()


class Main(QMainWindow, Ui_MainWindow):
    
    # Store miscellaneuous attributes here
    _nc = None
    _nc_groups = None
    _nc_group = None
    _nc_variables = None
    _nc_variable = None
    _cbar = None
    _use_cbar = False

    # Default options
    t = 0
    current_path = None
    timewindow = 1.00 # Plot window in days
    clevs = None
    cstep = 1. # Contour step
    cmap = 'Spectral_r'

    excluded_vars = ['zhat']

    def __init__(self,fig,ax ):
        super(Main, self).__init__()
        self.setupUi(self)

        # Create the figure and axes
        self.addmpl(fig,ax)
        # Add an open file menu bar option
        openFile = QAction(QtGui.QIcon('open.png'), 'Open NetCDF Station File', self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip('Open new File')
        openFile.triggered.connect(self.openfileDialog)

        menubar = self.menuBar()
        fileMenu = menubar.addMenu('&File')
        fileMenu.addAction(openFile)       

        ####
        self.init_cmapcombo()
        self.init_checkbox()
        self.init_windowtext()

        ###
        # 
        #self.ncfile = '/home/suntans/Share/ARCHub/DATA/FIELD/ShellPreludeRPS/NetCDF/Prelude_Gridded_Tuv.nc'
        ##self.ncfile = '/home/suntans/Share/ScottReef/DATA/FIELD/FALKOR/ProcessedData/FK150410_Gridded_Mooring_New.nc'

        #self.ncfile = '/home/suntans/Share/ScottReef/DATA/FIELD/FALKOR/ProcessedData/FK150410_Gridded_Mooring_New_TP.nc'

        # Open a file
        #self.open_ncfile()
        ## End testing
        ###
            
    def openfileDialog(self):

        if self.current_path is None:
            self.current_path = os.getcwd()

            fname, dstr = QFileDialog.getOpenFileName(self,\
            'Open station file', 
                    self.current_path,
            "NetCDF File (*.nc)")
            
        self.ncfile = str(fname)
        
        # Open a file
        self.open_ncfile()

    def open_ncfile(self):
        """
        Open the netcdf file and list all of the groups if any
        """
        if self._nc is not None:
            self._nc.close()

        try:
            self._nc = Dataset(self.ncfile, 'r')
        except:
            print('Failed to open file: %s'%self.ncfile)
            return

        # List the groups
        self._nc_groups = list(self._nc.groups.keys())

        # Print the groups in the Group ComboBox
        self.init_groupcombo()
        # Manually set to first group
        self.on_groupcombochange(0)
    
    def load_OceanMooring(self):
        """
        Load an ocean mooring object
        """
        print('Loading an OceanMooring object...')
        self.OM = oceanmooring.from_netcdf(\
            self.ncfile, self._nc_variable, group=self._nc_group)

        # Update the slider after the file is loaded
        self.init_slider()

        # Get the contours levels
        self.get_clevs()
    
    def init_groupcombo(self):
        cb = self.ncgroup_comboBox
        cb.clear()
        cb.addItems(self._nc_groups)
        cb.currentIndexChanged.connect(self.on_groupcombochange)
  		
    def init_variablecombo(self):
        cb = self.ncvariable_comboBox
        cb.clear()
        cb.addItems(self._nc_variables)
        cb.currentIndexChanged.connect(self.on_variablecombochange)

    def init_cmapcombo(self):
        cb = self.comboBox_cmap
        cb.clear()
        cb.addItems(cmaps)
        cb.currentIndexChanged.connect(self.on_cmapcombochange)

    def init_checkbox(self):
        b1 = self.checkBox_clim
        b1.setChecked(False)
        b1.stateChanged.connect(self.on_changecheckbox)

    def init_windowtext(self):
        # Do nothing
        b1 = self.textEdit_window
        #b1.triggered.connect(self.on_changewindow)
 
    def on_changecheckbox(self):
        b1 = self.checkBox_clim
        if b1.isChecked() == True:
            if self.timewindow != float(self.textEdit_window.toPlainText()):
                self.timewindow = float(self.textEdit_window.toPlainText())
            # Need to change the slider
            self.init_slider()

            # Re-draw the plot
            self.get_clevs()
            self.update_plot(self.t)
        #else:
        #    print 'Checkbox false'
    	
    def on_groupcombochange(self,i):
        self._nc_group =self._nc_groups[i]
        print("Selecting group: %s"%self._nc_group)
        self.list_variables(self._nc.groups[self._nc_group])

        # Update the first variable
        self.on_variablecombochange(0)
        #self._nc_variable = self._nc_variables[0]
        #self.load_OceanMooring()
        #self.update_plot(0)

    def on_variablecombochange(self,i):
        self._nc_variable =self._nc_variables[i]
        print("Selecting variable: %s"%self._nc_variable)

        # Now load the oceanmooring class
        self.load_OceanMooring()
        self.update_plot(0)

    def on_cmapcombochange(self,i):

        if USECMOCEAN:
            self.cmap = getattr(cm,cmaps[i])
        else:
            self.cmap = cmaps[i]

        self.update_plot(self.t)

    def on_timewindowchange(self):
        """

        """
        self.timewindow = float(self.textEdit_window.toPlainText())
        # Need to change the slider
        self.init_slider()

        # Update the plot
        self.get_clevs()
        self.update_plot(0)

    	
    def list_variables(self, ncobj):
        """
        List the variables in that group
        """
        self._nc_variables = []
        print(list(ncobj.variables.keys()))
        for vv in list(ncobj.variables.keys()):
            if hasattr(ncobj.variables[vv], 'coordinates') and \
                vv not in self.excluded_vars:
                self._nc_variables.append(vv)
        
        self.init_variablecombo()

    def init_slider(self):
        '''
        Update the slider based on the netcdf file
        '''
        #print dir(self.horizontalSlider)

        sl = self.horizontalSlider

        # Work out the time window size
        nt = self.OM.Nt

        windowsize = int(self.timewindow*86400.//self.OM.dt)
        self.windowsize = windowsize

        # Set these based on the time step
        sl.setMinimum(windowsize//2)
        sl.setMaximum(nt - windowsize - 1)
        sl.setValue(windowsize//2)
        sl.setTickPosition(QSlider.TicksBelow)
        sl.setTickInterval(windowsize)

        #layout.addWidget(self.sl)
        sl.valueChanged.connect(self.slider_valuechange)

    def slider_valuechange(self):
        self.t = self.horizontalSlider.value()
        self.update_plot(self.t)

    def keyPressEvent(self, event):
        """
        Change the behaviour of the left/right cursor
        """
        slider = self.horizontalSlider
        if event.key()==Qt.Key_Right:
            val = min(self.OM.Nt - self.windowsize//2,\
                slider.value() + self.windowsize//4)
            slider.setValue(val)
        elif event.key()==Qt.Key_Left:
            val = max(self.windowsize//2, slider.value() - self.windowsize//4)
            slider.setValue(val)
        else:
            QWidget.keyPressEvent(self, event)

    def update_plot(self, t):
        t0 = t
        t1 = t+self.windowsize

        if not self._cbar is None:
            self._use_cbar=False
        else:
            self._use_cbar=True
            
        # Update the contour plot
        self.ax.collections=[]
        C, cbar = self.OM.contourf(self.clevs, ax=self.ax,\
                fig=self.fig, t0=t0, t1=t1,
		filled=True, cbar=self._use_cbar, cmap=self.cmap)

        # update the colorbar
        if cbar is not None:
            self._cbar = cbar
        else:
            self._cbar.ax.cla()
            self._cbar = self.fig.colorbar(C, cax=self._cbar.ax)
            #print self._cbar, dir(self._cbar)
            #print dir(self.fig)
            #self._cbar.mappable = C
            #self._cbar.update_bruteforce(C)
            #self._cbar.on_mappable_changed(C)
            #self._cbar.update_normal(C)
            #self._cbar.set_clim([self.clevs[0], self.clevs[-1]])
            #self._cbar.update_ticks()
            #self._cbar.draw_all()


        # Add contours as well
        self.OM.contourf(self.clevs, ax=self.ax, t0=t0, t1=t1,
            filled=False, cbar=False,
            linecolor='0.5',linewidths=0.25)


        # Update the title
        #t0str = datetime.strftime(self.OM.t[t0], '%Y-%m-%d %H:%M:%S')
        #t1str = datetime.strftime(self.OM.t[t1], '%Y-%m-%d %H:%M:%S')
        t0str = np.datetime_as_string(self.OM.t[t0])[0:19].replace('T',' ')
        t1str = np.datetime_as_string(self.OM.t[t1])[0:19].replace('T',' ')

        self.ax.set_title('%s\n%s [%s]\n%s --> %s'%\
            (self.ncfile,self.OM.long_name, self.OM.units, t0str,t1str))

        #self.fig.tight_layout()
        self.fig.canvas.draw_idle()

    def get_clevs(self):
        """
        Change the contour levels here
        """
        if self.checkBox_clim.isChecked() == False:
            self.cmin = self.OM.y.min()
            self.cmax = self.OM.y.max()
            
            self.textEdit_cmin.setPlainText('%3.3f'%self.cmin)
            self.textEdit_cmax.setPlainText('%3.3f'%self.cmax)

        else:
            # Read the color limits from the data
            self.cmin = float(self.textEdit_cmin.toPlainText())
            self.cmax = float(self.textEdit_cmax.toPlainText())
            self.cstep = float(self.textEdit_contour.toPlainText())

            #
            self.cmax += self.cstep

        self.clevs = np.arange(self.cmin, self.cmax, self.cstep)

    def addmpl(self, fig, ax):
        """
        Add a matplotlib gui to the gui
        """
        self.fig = fig
        self.ax = ax
        self.canvas = FigureCanvas(fig)
        self.mplvl.addWidget(self.canvas)

        # Add the navigation toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.mplvl.addWidget(self.toolbar)

        self.canvas.draw()
		     
if __name__ == '__main__':
    import sys
    #from PyQt4.QtGui import QApplication
    from PyQt5.QtWidgets import QApplication
             
    # Create a matplotlib figure
    fig1 = Figure()
    ax1 = fig1.add_subplot(111)

    app = QApplication(sys.argv)
    main = Main(fig1, ax1)

    # Add the figure here
    #main.addmpl(fig1, ax1)
    main.show()
    sys.exit(app.exec_())
