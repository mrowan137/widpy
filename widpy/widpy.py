#!/usr/bin/env python
"""
widpy.py - Visualization tool for TRISTAN-MP data, with features specialized to 
simulations of magnetic reconnection.

Copyright 2018 Michael Rowan
Distributed under GPLv3 license.  See LICENSE for more information.

Based on the ImageView class, part of the pyqtgraph package, by Luke Campagnola.

Main features:
    *ROI mode for plotting spectra of particles
    *Rec mode to delineate reconnection region (requires fld files for right-tagged
     particles); threshold is controlled with slRec slider, which sets `dthresh`
    *Several matplotlib colorbars (viridis, RdBu, plasma, cubehelix)
"""
# Standard imports
import os
import sys
import time
import fnmatch
import numpy as np

# Matplotlib colormap
from matplotlib import cm

# PyQt5 Widgets
from PyQt5.QtWidgets import QApplication, QDialog
from PyQt5.QtWidgets import QFrame, QSplitter,  QVBoxLayout,  QHBoxLayout, QLabel
from PyQt5.QtWidgets import QSlider, QTextEdit, QLineEdit, QPushButton, QComboBox
from PyQt5.QtGui import QIcon
from PyQt5.QtCore import Qt

# Qt based plotting widgets
import pyqtgraph as pg

# For HDF5 files
import h5py

# Helper class, customized version of pyqtgraph ImageView
try:
    from helpers import ImageView
except ModuleNotFoundError:
    from widpy.helpers import ImageView


class Wid(QDialog):
    """
    Window containing three modified ImageViews, and information about load file.
    LineEdit allows for navigation to a different directory.
    """
    def __init__(self, parent=None, fileRoot=os.getcwd().replace('/output', '') + ('/')):
        super(Wid, self).__init__(parent)
        self.fileRoot = fileRoot
        self.initUI()

    def initUI(self):        
        """
        Load input data and Initialize widgets
        """
        # Read input, set constants, and read fluid file names
        self.getInput()
        self.setConstants()
        self.getKeys()
        
        # Initialize first ImageView
        # Set and pass size of the ROI `region of interest`
        roiSizey = self.paramDict['my0']/self.paramDict['istep']/7.
        roiSizex = roiSizey*self.golden
        self.imv1 = ImageView(roiPos=(0., 0.), roiSize=(roiSizex, roiSizey))        
        self.img1 = self.imv1.getImageItem()
        
        # Initialize a matplotlib colortable
        self.hist1 = self.imv1.getHistogramWidget()
        self.cbarName = 'viridis'
        self.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist1.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist1.gradient.ticks:
            tick.hide()
                                
        # Second ImageView
        self.imv2 = ImageView(roiPos=(0., 0.), roiSize=(roiSizex, roiSizey))                      
        self.hist2 = self.imv2.getHistogramWidget()
        self.cbarName = 'RdBu'
        self.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist2.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist2.gradient.ticks:
            tick.hide()

        # Third ImageView
        self.imv3 = ImageView(roiPos=(0., 0.), roiSize=(roiSizex, roiSizey))        
               
        self.hist3 = self.imv3.getHistogramWidget()
        self.cbarName = 'cubehelix'
        self.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist3.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist3.gradient.ticks:
            tick.hide()
            
        # Dropdown boxes to select fld file
        self.comboBox1 = QComboBox(self)
        self.comboBox2 = QComboBox(self)
        self.comboBox3 = QComboBox(self)
                
        # Fill dropdown boxes with fld file names
        for k in self.fldKeys:
            self.comboBox1.addItem(k)
            self.comboBox2.addItem(k)
            self.comboBox3.addItem(k)      

        self.comboBox1.activated[str].connect(self.onCombo1)        
        self.comboBox2.activated[str].connect(self.onCombo2)
        self.comboBox3.activated[str].connect(self.onCombo3)

        # Label, to list output step
        self.l1 = QLabel("          "+" ",self)
        self.l1.setAlignment=(Qt.AlignCenter)
        
        # Label, to list FPS
        self.l2 = QLabel("          "+" ",self)
        self.l2.setAlignment=(Qt.AlignCenter)
        
        # Label, to list dthresh
        self.l3 = QLabel("          "+" ",self)
        self.l3.setAlignment=(Qt.AlignCenter)
        
        # Label, to list prtl downsample factor
        self.l4 = QLabel("          "+" ",self)
        self.l4.setAlignment=(Qt.AlignCenter)
        
        # Slider to control file selection
        self.sl = QSlider(Qt.Horizontal)
        self.sl.setTracking(0)
        self.sl.setMinimum(int(self.fileRange[0]))
        self.sl.setMaximum(int(self.fileRange[-1]))
        self.sl.setValue(int(self.fileRange[0]))
        self.sl.setTickPosition(QSlider.TicksBelow)
        self.sl.setTickInterval(5)
        self.sl.valueChanged[int].connect(self.onActivateSlider)
        
        # Slider to control framerate
        self.slFPS = QSlider(Qt.Horizontal)
        self.slFPS.setMinimum(-12)
        self.slFPS.setMaximum(12)
        self.slFPS.setValue(1)
        self.slFPS.setTickPosition(QSlider.TicksBelow)
        self.slFPS.setTickInterval(1)
        self.slFPS.valueChanged[int].connect(self.onActivateSliderFPS)
        
        # Slider to set rec threshold
        self.slRec = QSlider(Qt.Horizontal)
        self.slRec.setTracking(0)
        self.slRec.setMinimum(1)
        self.slRec.setMaximum(50)
        self.slRec.setValue(30)
        self.slRec.setTickPosition(QSlider.TicksBelow)
        self.slRec.setTickInterval(1)
        self.slRec.valueChanged[int].connect(self.onActivateSliderRec)
        
        # Slider to downsample the prtl files
        self.slPrtlDown = QSlider(Qt.Horizontal)
        self.slPrtlDown.setTracking(0)
        self.slPrtlDown.setMinimum(1)
        self.slPrtlDown.setMaximum(20)
        self.slPrtlDown.setValue(10)
        self.slPrtlDown.setTickPosition(QSlider.TicksBelow)
        self.slPrtlDown.setTickInterval(1)
        self.slPrtlDown.valueChanged[int].connect(self.onActivateSliderPrtlDown)
        
        # Button, step left in time 
        self.TB1 = QPushButton("<", self)
        self.TB1.clicked.connect(self.onActivateTB1)
        
        # Button, step right in time
        self.TB2 = QPushButton(">", self)
        self.TB2.clicked.connect(self.onActivateTB2)
        
        # Interval for stepping through files
        self.tstep=1
        self.stepLine = QLineEdit(self)
        self.stepLine.setText("1")
        self.stepLine.returnPressed.connect(self.onReturn)
        
        # Navigation bar to select directory
        self.rootLine = QLineEdit(self)
        self.rootLine.setText(self.fileRoot)
        self.rootLine.returnPressed.connect(self.onReturnRoot)
        
        # List for paremeters (read from input file)
        self.paramList = QTextEdit(self)
        self.paramList.setReadOnly(True)
        for s in self.params:
            self.paramList.insertPlainText(s + '='+str(self.paramDict[s]))
            self.paramList.insertPlainText('\n')
        
        # Set the layout
        self.vbox = QVBoxLayout(self) 
        self.space = QFrame(self)
        self.space.setFrameShape(QFrame.StyledPanel)
        
        # Splitters, to allow for resizable frames
        self.splitterMV = QSplitter(Qt.Vertical)
        self.splitterMH = QSplitter(Qt.Horizontal)
        self.splitterL = QSplitter(Qt.Vertical)        
        self.splitterR = QSplitter(Qt.Vertical)
        self.splitterTButtons = QSplitter(Qt.Horizontal)
        
        # ImageViews go on the right
        self.splitterR.addWidget(self.imv1)
        self.imv1.parent = self
        self.splitterR.addWidget(self.imv2)
        self.imv2.parent = self
        self.splitterR.addWidget(self.imv3)
        self.imv3.parent = self
        
        # Selection panels go on the left
        self.splitterL.addWidget(self.comboBox1)
        self.splitterL.addWidget(self.comboBox2)
        self.splitterL.addWidget(self.comboBox3)

        # Label and sliders on the left
        self.splitterL.addWidget(self.l1)
        self.splitterL.addWidget(self.sl)
        
        # Navigation buttons on the left
        self.splitterTButtons.addWidget(self.TB1)
        self.splitterTButtons.addWidget(self.stepLine)
        self.splitterTButtons.addWidget(self.TB2)
        self.splitterTButtons.setSizes([10, 5, 10])
        self.splitterL.addWidget(self.splitterTButtons)
        
        # FPS label and slider on the left
        self.splitterL.addWidget(self.l2)
        self.splitterL.addWidget(self.slFPS)

        # Rec label and slider on the left
        self.splitterL.addWidget(self.l3)
        self.splitterL.addWidget(self.slRec)    
        
        # Rec label and slider on the left
        self.splitterL.addWidget(self.l4)
        self.splitterL.addWidget(self.slPrtlDown)   
    
        # Parameter list on the left
        self.splitterL.addWidget(self.paramList)
        self.splitterL.addWidget(self.space)
        self.splitterL.setSizes([10, 10, 10, 10, 10, 10, 10,
                                 10, 10, 10, 10, 10, 730, 10])
        
        # Split ImageViews and time navigation tools 
        self.splitterMH.addWidget(self.splitterL)
        self.splitterMH.addWidget(self.splitterR)
        self.splitterMH.setSizes([250, 550])
        
        # Split directory navigation from rest
        self.splitterMV.addWidget(self.rootLine)
        self.splitterMV.addWidget(self.splitterMH)
        self.splitterMV.setSizes([30, 770])
        
        self.vbox.addWidget(self.splitterMV)
        self.setLayout(self.vbox)
        self.setGeometry(1000, 0, 800, 800)
        self.setWindowTitle('widpy')
        
        # Activate panels
        self.l1.setText('output file: ' + str(self.sl.value()) 
                        + '/' + str(self.fileRange[-1]))

        self.fps = self.slFPS.value()
        self.l2.setText('FPS: ' + str(self.slFPS.value()))

        self.slRec.setEnabled(False)
        self.dthresh = 0.01*self.slRec.value()
        self.l3.setText('dthresh: ' + str(self.dthresh))
        
        self.prtlDown = self.slPrtlDown.value()
        self.l4.setText('Prtl downsample factor: ' + str(self.prtlDown)) 

        self.text1 = 'dens'
        self.text2 = 'v3y'
        self.text3 = 'jz'
    
        self.comboBox1.setCurrentText(self.text1)
        self.comboBox2.setCurrentText(self.text2)
        self.comboBox3.setCurrentText(self.text3)
        
        # zrange for plots, colorbar        
        self.onCombo1(self.text1)
        self.onCombo2(self.text2)
        self.onCombo3(self.text3)
           
    def setCbar(self):
        """
        Make colorbars from matplotlib
        """
        # Get the colormap, select rgb values
        colormap = cm.get_cmap(self.cbarName) 
        colormap._init()
        lut = (colormap._lut*1).view(np.ndarray)
        self.colorTable = np.array(lut[:,0:-1])
            
    def getKeys(self):
        """
        Get the names of fluid files
        """
        self.fileload = self.fileRoot+self.fileEndFld
        f = h5py.File(self.fileload, 'r')
        self.fldKeys = []
        for k in f.keys():
            self.fldKeys = self.fldKeys + [k]
            
        # Add electron keys
        self.fldKeysLecs = ['dense', 'v3xe', 'v3ye', 'v3ze']
        for k in self.fldKeysLecs:
            self.fldKeys = self.fldKeys + [k]
            
        # Add energy keys
        self.fldKeysUint = ['keneint', 'keniint', 'bulke', 'bulki']
        for k in self.fldKeysUint:
            self.fldKeys = self.fldKeys + [k]
            
        # Add E x B drift velocity keys
        self.fldKeysDrift = ['vdrx', 'vdry', 'vdrz']
        for k in self.fldKeysDrift:
            self.fldKeys = self.fldKeys + [k]
            
        # Add curl B keys
        self.fldKeysCurlB = ['crlbx', 'crlby', 'crlbz']
        for k in self.fldKeysCurlB:
            self.fldKeys = self.fldKeys + [k]

        # Add curl E keys
        self.fldKeysCurlE = ['crlex', 'crley', 'crlez']
        for k in self.fldKeysCurlE:
            self.fldKeys = self.fldKeys + [k]
            
        # Add jdotE key
        self.fldKeys = self.fldKeys + ['jdotE'] 
        
        # Add non-ideal E field, elec
        self.fldKeysNonidealElec = ['ex_nonideal_e', 'ey_nonideal_e', 'ez_nonideal_e']
        for k in self.fldKeysNonidealElec:
            self.fldKeys = self.fldKeys + [k]

        # Add non-ideal E field, ions
        self.fldKeysNonidealIon = ['ex_nonideal_i', 'ey_nonideal_i', 'ez_nonideal_i']
        for k in self.fldKeysNonidealIon:
            self.fldKeys = self.fldKeys + [k]
        
    def getInput(self):
        """
        Use current directory as load file (else use default), and read input
        """
        # check for output directory; use default data if not found
        if os.path.exists(self.fileRoot+'/output'):
            dr = os.listdir(self.fileRoot+'/output')
        else:
            print('No TRISTAN output directory found!  Using default data.') 
            self.fileRoot = (os.path.dirname(os.path.abspath(__file__)) 
                  + '/../tests/test_data/')
            dr = os.listdir(self.fileRoot+'/output')
        
        # Get the min/max of output files
        self.fileRange = np.array([])
        for f in dr:
            if fnmatch.fnmatch(f, 'flds.tot.*'):
                self.fileRange = np.append(self.fileRange, f[-3::])
            
        if len(self.fileRange) > 0:
            # Sort the filelist
            argsrt = np.argsort(self.fileRange)
            self.fileRange = self.fileRange[argsrt]
            
            # Set start, stop
            self.fileStartFld = 'output/flds.tot.' + self.fileRange[0]
            self.fileStartPrtl = 'output/prtl.tot.' + self.fileRange[0]
            self.fileEndFld = 'output/flds.tot.' + self.fileRange[-1]
            self.fileEndPrtl = 'output/prtl.tot.' + self.fileRange[-1]
                       
        # Read the input file
        with open(self.fileRoot + "input", 'r') as f:
            contents = f.readlines()
            
        params = []
        vals = []
        str1, str2, str3 = ' '*3

        for line in contents:
            for s in line.split():
                str1 = str2
                str2 = str3
                str3 = s
                if (str2 == "=" and str1[0] != "#"): 
                    params.append(str1)
                    vals.append(float(str3))

        z = zip(params,vals)
        paramDict = dict(z)
        self.paramDict = paramDict
        self.z = z
        self.params = params
        self.vals = vals
        f.close()
        
    def setConstants(self):
        """
        Set useful constants
        """
        # Used in plotting
        self.golden = self.golden = (1 + 5**0.5)/2
                        
    def onReturnRoot(self):        
        """
        Get root entered to navigation bar, reinitialize
        """
        text = self.rootLine.text()
        self.fileRoot = text
        self.rootLine.setText(text)
        self.fileLoad = text + self.fileEndFld[0:-3] + self.strn1z
        print('Changing directory: ',self.fileRoot)
        
        # Get fld keys, input, update layout
        self.getInput()
        self.getKeys()
        self.updateLayout()

        # Hack
        self.sl.setValue(self.sl.value()+(self.tstep))
        
    def onReturn(self):
        """
        Set interval to step through files when return key is pressed (for LineEdit
        between L/R pressbuttons)
        """
        text = self.stepLine.text()
        self.tstep = int(text)
        self.sl.setValue(self.sl.value() + (self.tstep))

    def onActivated1(self, text):
        """
        Load data, update ImageView1
        """
        # Load data, flds
        self.text1 = text
        self.strn1 = str(self.sl.value())
        self.strn1z = self.strn1.zfill(3)

        # If electron quantity, need to compute from ion and total
        if (text == 'dense'):
            # Electron number density
            datasetFull = (self.loadDataFld('dens', self.strn1z) - 
                           self.loadDataFld('densi', self.strn1z))
        elif (text in ['v3xe', 'v3ye', 'v3ze']):
            # x, y, z components of electron fluid velocity/c
            v3j = self.loadDataFld(text.replace('e',''), self.strn1z)
            v3ji = self.loadDataFld(text.replace('e','i'), self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            dense = dens - densi
            v3je = ((v3j*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3ji*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))

            # Electron 3-velocity is zero if no particles
            v3je[dense == 0.] = 0.
            datasetFull = v3je
        elif (text == 'keneint'):
            # Electron dimensionless internal energy u_int/(m_e c^2);
            # Calculation here is appropriate only in non-relativistic limit
            ken = self.loadDataFld('ken', self.strn1z)
            keni = self.loadDataFld('keni', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            kene = ((ken*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - keni*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            keneint = kene - 1 - bulke
            datasetFull = keneint
        elif (text == 'keniint'):
            # Ion dimensionless internal energy u_int/(m_i c^2);
            # Calculation here is appropriate only in non-relativistic limit
            keni = self.loadDataFld('keni', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)            
            keniint = keni - 1 - bulki
            datasetFull = keniint
        elif (text == 'bulki'):
            # Ion dimensionless bulk energy 0.5*(v_i/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            datasetFull = bulki
        elif (text == 'bulke'):
            # Electron dimensionless bulk energy 0.5*(v_e/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            datasetFull = bulke
        elif (text == 'vdrx'):
            # E x B drift velocity, x component
            ey = self.loadDataFld('ey', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, x component
            vdrx = (ey*bz - ez*by)/bsq

            # Cells where B = 0
            vdrx[bsq == 0.] = 0.
            datasetFull = vdrx
        elif (text == 'vdry'):
            # E x B drift velocity, y component
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            ex = self.loadDataFld('ex', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, y component
            vdry = (ez*bx - ex*bz)/bsq

            # Cells where B = 0
            vdry[bsq == 0.] = 0.
            datasetFull = vdry
        elif (text == 'vdrz'):
            # E x B drift velocity, z component
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, z component
            vdrz = (ex*by - ey*bx)/bsq

            # Cells where B = 0
            vdrz[bsq == 0.] = 0.
            datasetFull = vdrz
        elif (text == 'crlbx'):
            # Curl of B, x component
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            dby_dz = 0.
            dbz_dy, _ = np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbx = dby_dz - dbz_dy
            datasetFull = crlbx
        elif (text == 'crlby'):
            # Curl of B, y component
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            _, dbz_dx= np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dbx_dz = 0.

            crlby = dbz_dx - dbx_dz
            datasetFull = crlby
        elif (text == 'crlbz'):
            # Curl of B, z component
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            dbx_dy, _ = np.gradient(bx, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dby_dx= np.gradient(by, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbz = dbx_dy - dby_dx
            datasetFull = crlbz
        elif (text == 'crlex'):
            # Curl of E, x component
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            dey_dz = 0.
            dez_dy, _ = np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])

            crlex = dey_dz - dez_dy
            datasetFull = crlex
        elif (text == 'crley'):
            # Curl of E, y component
            ex = self.loadDataFld('ex', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            _, dez_dx= np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dex_dz = 0.

            crley = dez_dx - dex_dz
            datasetFull = crley
        elif (text == 'crlez'):
            # Curl of E, z component
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            dex_dy, _ = np.gradient(ex, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dey_dx= np.gradient(ey, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])

            crlez = dex_dy - dey_dx
            datasetFull = crlez
        elif (text == 'jdotE'):
            # Energization j.E, power/unit volume
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            jx = self.loadDataFld('jx', self.strn1z)
            jy = self.loadDataFld('jy', self.strn1z)
            jz = self.loadDataFld('jz', self.strn1z)
            
            jdotE = ex*jx + ey*jy + ez*jz
            datasetFull = jdotE
        elif (text == 'ex_nonideal_e'):
            # Non-ideal E field, x component, electrons
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ex_nonideal_e = ex + (v3ye*bz - v3ze*by)
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_e'):
            # Non-ideal E field, y component, electrons
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ey_nonideal_e = ey + (v3ze*bx - v3xe*bz)
            datasetFull = ey_nonideal_e
        elif (text == 'ez_nonideal_e'):
            # Non-ideal E field, z component, electrons
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.

            ez_nonideal_e = ez + (v3xe*by - v3ye*bx)
            datasetFull = ez_nonideal_e
        elif (text == 'ex_nonideal_i'):
            # Non-ideal E field, z component, ions
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ex_nonideal_i = ex + (v3yi*bz - v3zi*by)            
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_i'):
            # Non-ideal E field, y component, ions
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ey_nonideal_i = ey + (v3zi*bx - v3xi*bz)
            datasetFull = ey_nonideal_i
        elif (text == 'ez_nonideal_i'):
            # Non-ideal E field, z component, ions
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)

            ez_nonideal_i = ez + (v3xi*by - v3yi*bx)
            datasetFull = ez_nonideal_i
        else:
            datasetFull = self.loadDataFld(text, self.strn1z)

        dataset = self.subsetData(datasetFull)

        # Load data, prtls if available; convert the x-coordinate
        if self.imv1.roiMode:
            self.xe = self.loadDataPrtl('xe', self.strn1z)/self.paramDict['istep'] + self.xshift
            self.ye = self.loadDataPrtl('ye', self.strn1z)/self.paramDict['istep']
            self.gammae = self.loadDataPrtl('gammae', self.strn1z)
        
            # Pass prtl data to the imageview
            self.imv1.setScatter((self.xe, self.ye, self.gammae))
        
        if self.imv1.recMode:
            pdensrFull = self.loadDataFld('pdens', self.strn1z)
            pdensr = self.subsetData(pdensrFull)
            densFull = self.loadDataFld('dens', self.strn1z)
            dens = self.subsetData(densFull)
            bdensFull = self.loadDataFld('bdens', self.strn1z)
            bdens = self.subsetData(bdensFull) 
            
            pdenstot = dens - bdens
            
            rat = self.nandiv(pdensr, pdenstot)
            whrec=np.where((rat > self.dthresh) & 
                                (rat < (1-self.dthresh)) &
                                (bdens == 0))
            
            my, mx = np.shape(bdens)
            self.recReg = np.zeros((my, mx))
            self.recReg[whrec]=1.
            
        # Update image
        levs = self.imv1.getLevels()
        self.imv1.setImage(dataset, autoHistogramRange=False, autoRange=False)
        self.imv1.setLevels(min=levs[0], max=levs[-1])
        self.hist1.setHistogramRange(levs[0], levs[-1], padding=0.1)
        
    def onActivated2(self,text):        
        """
        Load data, update ImageView2
        """
        # Load data, flds
        self.text2 = text
        self.strn2 = str(self.sl.value())
        self.strn2z = self.strn2.zfill(3)

        # If electron quantity, need to compute from ion and total
        if (text == 'dense'):
            # Electron number density
            datasetFull = (self.loadDataFld('dens', self.strn2z) - 
                           self.loadDataFld('densi', self.strn2z))
        elif (text in ['v3xe', 'v3ye', 'v3ze']):
            # x, y, z components of electron fluid velocity/c
            v3j = self.loadDataFld(text.replace('e',''), self.strn2z)
            v3ji = self.loadDataFld(text.replace('e','i'), self.strn2z)
            dens = self.loadDataFld('dens', self.strn2z)
            densi = self.loadDataFld('densi', self.strn2z)
            dense = dens - densi
            v3je = ((v3j*(densi+(1./(self.paramDict['mi']/self.paramDict['me']))*dense)
                     -v3ji*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))

            # Electron 3-velocity is zero if no particles
            v3je[dense == 0.] = 0.            
            datasetFull = v3je
        elif (text == 'keneint'):
            # Electron dimensionless internal energy u_int/(m_e c^2);
            # Calculation here is appropriate only in non-relativistic limit
            ken = self.loadDataFld('ken', self.strn1z)
            keni = self.loadDataFld('keni', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            kene = ((ken*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense)
                     - keni*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)            
            keneint = kene - 1 - bulke
            datasetFull = keneint
        elif (text == 'keniint'):
            # Ion dimensionless internal energy u_int/(m_i c^2);
            # Calculation here is appropriate only in non-relativistic limit
            keni = self.loadDataFld('keni', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            keniint = keni - 1 - bulki
            datasetFull = keniint
        elif (text == 'bulki'):
            # Ion dimensionless bulk energy 0.5*(v_i/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            datasetFull = bulki
        elif (text == 'bulke'):
            # Electron dimensionless bulk energy 0.5*(v_e/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            datasetFull = bulke
        elif (text == 'vdrx'):
            # E x B drift velocity, x component
            ey = self.loadDataFld('ey', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, x component
            vdrx = (ey*bz - ez*by)/bsq

            # Cells where B = 0
            vdrx[bsq == 0.] = 0.            
            datasetFull = vdrx
        elif (text == 'vdry'):
            # E x B drift velocity, y component
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            ex = self.loadDataFld('ex', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, y component
            vdry = (ez*bx - ex*bz)/bsq

            # Cells where B = 0
            vdry[bsq == 0.] = 0.
            datasetFull = vdry
        elif (text == 'vdrz'):
            # E x B drift velocity, z component
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, z component
            vdrz = (ex*by - ey*bx)/bsq

            # Cells where B = 0
            vdrz[bsq == 0.] = 0.
            datasetFull = vdrz
        elif (text == 'crlbx'):
            # Curl of B, x component
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            dby_dz = 0.
            dbz_dy, _ = np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbx = dby_dz - dbz_dy            
            datasetFull = crlbx
        elif (text == 'crlby'):
            # Curl of B, y component
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            _, dbz_dx= np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dbx_dz = 0.

            crlby = dbz_dx - dbx_dz
            datasetFull = crlby
        elif (text == 'crlbz'):
            # Curl of B, z component
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            dbx_dy, _ = np.gradient(bx, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dby_dx= np.gradient(by, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbz = dbx_dy - dby_dx
            datasetFull = crlbz
        elif (text == 'crlex'):
            # Curl of E, x component
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            dey_dz = 0.
            dez_dy, _ = np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])

            crlex = dey_dz - dez_dy
            datasetFull = crlex
        elif (text == 'crley'):
            # Curl of E, y component
            ex = self.loadDataFld('ex', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            _, dez_dx= np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dex_dz = 0.

            crley = dez_dx - dex_dz
            datasetFull = crley
        elif (text == 'crlez'):
            # Curl of E, z component
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            dex_dy, _ = np.gradient(ex, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dey_dx= np.gradient(ey, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])

            crlez = dex_dy - dey_dx
            datasetFull = crlez
        elif (text == 'jdotE'):
            # Energization j.E, power/unit volume
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            jx = self.loadDataFld('jx', self.strn1z)
            jy = self.loadDataFld('jy', self.strn1z)
            jz = self.loadDataFld('jz', self.strn1z)
            
            jdotE = ex*jx + ey*jy + ez*jz
            datasetFull = jdotE
        elif (text == 'ex_nonideal_e'):
            # Non-ideal E field, x component, electrons
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ex_nonideal_e = ex + (v3ye*bz - v3ze*by)
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_e'):
            # Non-ideal E field, y component, electrons
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ey_nonideal_e = ey + (v3ze*bx - v3xe*bz)
            datasetFull = ey_nonideal_e
        elif (text == 'ez_nonideal_e'):
            # Non-ideal E field, z component, electrons
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.

            ez_nonideal_e = ez + (v3xe*by - v3ye*bx)
            datasetFull = ez_nonideal_e
        elif (text == 'ex_nonideal_i'):
            # Non-ideal E field, x component, ions
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ex_nonideal_i = ex + (v3yi*bz - v3zi*by)
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_i'):
            # Non-ideal E field, y component, ions
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ey_nonideal_i = ey + (v3zi*bx - v3xi*bz)
            datasetFull = ey_nonideal_i
        elif (text == 'ez_nonideal_i'):
            # Non-ideal E field, z component, ions
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            
            ez_nonideal_i = ez + (v3xi*by - v3yi*bx)
            datasetFull = ez_nonideal_i
        else:
            datasetFull = self.loadDataFld(text, self.strn2z)

        dataset = self.subsetData(datasetFull)
        
        # Load data, prtls if available; convert the x-coordinate   
        if self.imv2.roiMode:
            self.xe = self.loadDataPrtl('xe', self.strn2z)/self.paramDict['istep'] + self.xshift
            self.ye = self.loadDataPrtl('ye', self.strn2z)/self.paramDict['istep']
            self.gammae = self.loadDataPrtl('gammae', self.strn2z)
        
            # Pass prtl data to the imageview
            self.imv2.setScatter((self.xe, self.ye, self.gammae))
        
        # Update image
        levs = self.imv2.getLevels()
        self.imv2.setImage(dataset, autoHistogramRange=False, autoRange=False)
        self.imv2.setLevels(min=levs[0], max=levs[1])        
        self.hist2.setHistogramRange(levs[0], levs[-1], padding=0.1)

    def onActivated3(self,text):
        """
        Load data, update ImageView1
        """
        # Load data, flds
        self.text3 = text
        self.strn3 = str(self.sl.value())
        self.strn3z = self.strn3.zfill(3)

        # If electron quantity, need to compute from ion and total
        if (text == 'dense'):
            # Electron number density
            datasetFull = (self.loadDataFld('dens', self.strn3z) - 
                           self.loadDataFld('densi', self.strn3z))
        elif (text in ['v3xe', 'v3ye', 'v3ze']):
            # x, y, z components of electron fluid velocity/c
            v3j = self.loadDataFld(text.replace('e',''), self.strn3z)
            v3ji = self.loadDataFld(text.replace('e','i'), self.strn3z)
            dens = self.loadDataFld('dens', self.strn3z)
            densi = self.loadDataFld('densi', self.strn3z)
            dense = dens - densi
            v3je = ((v3j*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3ji*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))

            # Electron 3-velocity is zero if no particles
            v3je[dense == 0.] = 0.
            datasetFull = v3je
        elif (text == 'keneint'):
            # Electron dimensionless internal energy u_int/(m_e c^2);
            # Calculation here is appropriate only in non-relativistic limit
            ken = self.loadDataFld('ken', self.strn1z)
            keni = self.loadDataFld('keni', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            kene = ((ken*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - keni*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            keneint = kene - 1 - bulke
            datasetFull = keneint
        elif (text == 'keniint'):
            # Ion dimensionless internal energy u_int/(m_i c^2);
            # Calculation here is appropriate only in non-relativistic limit
            keni = self.loadDataFld('keni', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            keniint = keni - 1 - bulki
            datasetFull = keniint
        elif (text == 'bulki'):
            # Ion dimensionless bulk energy 0.5*(v_i/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            datasetFull = bulki
        elif (text == 'bulke'):
            # Electron dimensionless bulk energy 0.5*(v_e/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            datasetFull = bulke
        elif (text == 'vdrx'):
            # E x B drift velocity, x component
            ey = self.loadDataFld('ey', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, x component
            vdrx = (ey*bz - ez*by)/bsq

            # Cells where B = 0
            vdrx[bsq == 0.] = 0.            
            datasetFull = vdrx
        elif (text == 'vdry'):
            # E x B drift velocity, y component
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            ex = self.loadDataFld('ex', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, y component
            vdry = (ez*bx - ex*bz)/bsq

            # Cells where B = 0
            vdry[bsq == 0.] = 0.
            datasetFull = vdry
        elif (text == 'vdrz'):
            # E x B drift velocity, z component
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, z component
            vdrz = (ex*by - ey*bx)/bsq

            # Cells where B = 0
            vdrz[bsq == 0.] = 0.
            datasetFull = vdrz
        elif (text == 'crlbx'):
            # Curl of B, x component
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            dby_dz = 0.
            dbz_dy, _ = np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbx = dby_dz - dbz_dy            
            datasetFull = crlbx
        elif (text == 'crlby'):
            # Curl of B, y component
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            _, dbz_dx= np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dbx_dz = 0.
            
            crlby = dbz_dx - dbx_dz            
            datasetFull = crlby
        elif (text == 'crlbz'):
            # Curl of B, z component
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            dbx_dy, _ = np.gradient(bx, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dby_dx= np.gradient(by, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbz = dbx_dy - dby_dx
            datasetFull = crlbz
        elif (text == 'crlex'):
            # Curl of E, x component
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            dey_dz = 0.
            dez_dy, _ = np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlex = dey_dz - dez_dy
            datasetFull = crlex
        elif (text == 'crley'):
            # Curl of E, y component
            ex = self.loadDataFld('ex', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            _, dez_dx= np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dex_dz = 0.

            crley = dez_dx - dex_dz
            datasetFull = crley
        elif (text == 'crlez'):
            # Curl of E, z component
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            dex_dy, _ = np.gradient(ex, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dey_dx= np.gradient(ey, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlez = dex_dy - dey_dx
            datasetFull = crlez
        elif (text == 'jdotE'):
            # Energization j.E, power/unit volume
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            jx = self.loadDataFld('jx', self.strn1z)
            jy = self.loadDataFld('jy', self.strn1z)
            jz = self.loadDataFld('jz', self.strn1z)
            
            jdotE = ex*jx + ey*jy + ez*jz
            datasetFull = jdotE
        elif (text == 'ex_nonideal_e'):
            # Non-ideal E field, x component, electrons
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ex_nonideal_e = ex + (v3ye*bz - v3ze*by)
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_e'):
            # Non-ideal E field, y component, electrons
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ey_nonideal_e = ey + (v3ze*bx - v3xe*bz)
            datasetFull = ey_nonideal_e
        elif (text == 'ez_nonideal_e'):
            # Non-ideal E field, z component, electrons
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.

            ez_nonideal_e = ez + (v3xe*by - v3ye*bx)
            datasetFull = ez_nonideal_e
        elif (text == 'ex_nonideal_i'):
            # Non-ideal E field, x component, ions
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ex_nonideal_i = ex + (v3yi*bz - v3zi*by)
            datasetFull = ex_nonideal_e
        elif (text == 'ey_nonideal_i'):
            # Non-ideal E field, y component, ions
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ey_nonideal_i = ey + (v3zi*bx - v3xi*bz)
            datasetFull = ey_nonideal_i
        elif (text == 'ez_nonideal_i'):
            # Non-ideal E field, z component, ions
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)

            ez_nonideal_i = ez + (v3xi*by - v3yi*bx)
            datasetFull = ez_nonideal_i
        else:
            datasetFull = self.loadDataFld(text, self.strn3z)

        dataset = self.subsetData(datasetFull)
        # Load data, prtls if available; convert the x-coordinate        
        if self.imv3.roiMode:
            self.xe = self.loadDataPrtl('xe', self.strn3z)/self.paramDict['istep'] + self.xshift
            self.ye = self.loadDataPrtl('ye', self.strn3z)/self.paramDict['istep']
            self.gammae = self.loadDataPrtl('gammae', self.strn3z)
        
            # Pass prtl data to the imageview
            self.imv3.setScatter((self.xe, self.ye, self.gammae))
        
        # Update image
        levs = self.imv3.getLevels()
        self.imv3.setImage(dataset, autoHistogramRange=False, autoRange=False)
        self.imv3.setLevels(min=levs[0], max=levs[1]) 
        self.hist3.setHistogramRange(levs[0], levs[-1], padding=0.1)
           
    def onCombo1(self, t):
        """
        Same as onActivate1, but also compute the zrange
        """
        self.levsMin1, self.levsMax1 = self.getLevsMinMax(t)
        self.onActivated1(t)
        self.imv1.setLevels(min=self.levsMin1, max=self.levsMax1)
        self.hist1.setHistogramRange(self.levsMin1, self.levsMax1, padding=0.1)

    def onCombo2(self, t):
        """
        Same as onActivate2, but also compute the zrange
        """
        self.levsMin2, self.levsMax2 = self.getLevsMinMax(t)
        self.onActivated2(t)
        self.imv2.setLevels(min=self.levsMin2, max=self.levsMax2)
        self.hist2.setHistogramRange(self.levsMin2, self.levsMax2, padding=0.1)

    def onCombo3(self, t):
        """
        Same as onActivate3, but also compute the zrange
        """
        self.levsMin3, self.levsMax3 = self.getLevsMinMax(t)
        self.onActivated3(t)
        self.imv3.setLevels(min=self.levsMin3, max=self.levsMax3)
        self.hist3.setHistogramRange(self.levsMin3, self.levsMax3, padding=0.1)
        
    def onActivateSlider(self, t):
        """
        Set the current file index
        """
        self.currInd = np.where(t == np.array([int(j) for j in self.fileRange]))[0]
        
        self.l1.setText("output file: " + str(t) + '/' + str(self.fileRange[-1]))
        self.onActivated1(self.text1)
        self.onActivated2(self.text2)
        self.onActivated3(self.text3)
        
    def onActivateSliderFPS(self):
        """
        Set slider value to framerate (zero is not allowed)
        """
        if (self.slFPS.value() !=0):
            self.fps = self.slFPS.value()
            self.l2.setText("FPS: " + str(self.fps))
        elif (self.slFPS.value() == 0 and self.fps > 0):
            self.slFPS.setValue(-1)
            self.fps = -1
            self.l2.setText("FPS: " + str(-1))
        elif (self.slFPS.value() == 0 and self.fps < 0):
            self.slFPS.setValue(1)
            self.fps = 1
            self.l2.setText("FPS: " + str(1))
            
    def onActivateSliderRec(self):
        """
        Set slider value to reconnection threshold
        """
        self.dthresh = 0.01*self.slRec.value()
        self.l3.setText("dthresh: " + str(self.dthresh))
        self.onActivated1(self.text1)
        self.onActivated2(self.text2)
        self.onActivated3(self.text3)
        
    def onActivateSliderPrtlDown(self):
        """
        Set slider value to prtl file downsample factor
        """
        self.prtlDown = self.slPrtlDown.value()
        self.l4.setText("prtl downsample factor: " + str(self.prtlDown))
        if (self.imv1.roiMode):
            self.onActivated1(self.text1)
        if (self.imv2.roiMode):
            self.onActivated2(self.text2)
        if (self.imv3.roiMode):
            self.onActivated3(self.text3)
        
        
    def onActivateTB1(self):
        """
        Activation of left pressbutton to step through files
        """
        self.sl.setValue((self.sl.value() - self.tstep 
                          - int(self.fileRange[0]))%len(self.fileRange)
                         + int(self.fileRange[0]))

    def onActivateTB2(self):
        """
        Activation of right pressbutton to step through files
        """
        self.sl.setValue((self.sl.value() + self.tstep 
                          - int(self.fileRange[0]))%len(self.fileRange)
                         + int(self.fileRange[0]))

    def updateLayout(self):
        """
        Clear comboBoxes and update, along with slider and parameter list,
        according to the current fileroot
        """
        self.comboBox1.clear()
        self.comboBox2.clear()
        self.comboBox3.clear()
        for k in self.fldKeys:
            self.comboBox1.addItem(k)
            self.comboBox2.addItem(k)
            self.comboBox3.addItem(k)   
            
        self.text1 = self.comboBox1.itemText(0)
        self.text2 = self.comboBox2.itemText(1)
        self.text3 = self.comboBox3.itemText(2)
    
        self.comboBox1.setCurrentText(self.text1)
        self.comboBox2.setCurrentText(self.text2)
        self.comboBox3.setCurrentText(self.text3)

        # Slider
        self.sl.setMinimum(self.fileRange[0].astype(int))
        self.sl.setMaximum(self.fileRange[-1].astype(int))
        
        # Parameter list
        self.paramList.clear()
        for s in self.params:
            self.paramList.insertPlainText(s + '='+str(self.paramDict[s]))
            self.paramList.insertPlainText('\n')
                
    def loadDataFld(self, field, num):
        """
        Load an hdf5 file of format `fld.num`,  using the current self.fileRoot
           *field - (string) field to load, e.g. 'v3x', 'ken', ...
           *num - (string) 3 characters is length: '001', '002', ...
        """
        fileLoad = self.fileRoot + 'output/flds.tot.' + num
        f = h5py.File(fileLoad, 'r')
        data = f[field][0,:,:]
        return data
    
    def loadDataPrtl(self, field, num):
        """
        Load an hdf5 file of format `fld.num`,  using the current self.fileRoot
           *field - (string) field to load, e.g. 'xe', 'gammae', ...
           *num - (string) 3 characters is length: '001', '002', ...
        """
        fileLoad = self.fileRoot + 'output/prtl.tot.' + num
        f = h5py.File(fileLoad, 'r')
        data = f[field][::self.prtlDown]
        return data
        
    def subsetData(self, data):
        """
        Select a central region of a 2D array
        """
        my, mx = data.shape
        xmid = int(mx/2.)
        golden = (1 + 5 ** 0.5) / 2
        deltax = min(int(0.5*my/golden), int(0.5*mx))
        deltax = xmid
        subset = data[:,xmid - deltax:xmid + deltax]
        # Based on the selection, a shift to align prtl coordinate with the grid
        self.xshift = deltax - xmid
        
        return subset

    def nandiv(self, m1, m2):
        """
        Divid array m1 by b2, set resultant array to zero in cells where m2 is zero
        """
        tmp1 = m1
        whz = np.where(m2 == 0.)
        tmp2 = m2
        tmp2[whz]=1.
        res = tmp1/tmp2
        res[whz] = 0.
        return res
    
    def getLevsMinMax(self, text):
        """
        Get min/max levels for histogram and image
        """
        # If electron quantity, need to compute from ion and total
        if (text == 'dense'):
            # Electron number density
            data = (self.loadDataFld('dens', self.fileRange[-1]) - 
                    self.loadDataFld('densi', self.fileRange[-1]))
        elif (text in ['v3xe', 'v3ye', 'v3ze']):
            # x, y, z components of electron fluid velocity/c
            v3j = self.loadDataFld(text.replace('e',''), self.strn1z)
            v3ji = self.loadDataFld(text.replace('e','i'), self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            dense = dens - densi
            v3je = ((v3j*(densi+(1./(self.paramDict['mi']/self.paramDict['me']))*dense)-v3ji*densi)
                         /((1./(self.paramDict['mi']/self.paramDict['me']))*dense))

            # Electron 3-velocity is zero if no particles
            v3je[dense == 0.] = 0.
            data = v3je
        elif (text == 'keneint'):
            # Electron dimensionless internal energy u_int/(m_e c^2);
            # Calculation here is appropriate only in non-relativistic limit
            ken = self.loadDataFld('ken', self.strn1z)
            keni = self.loadDataFld('keni', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            kene = ((ken*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense)
                     - keni*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            keneint = kene - 1 - bulke
            data = keneint
        elif (text == 'keniint'):
            # Ion dimensionless internal energy u_int/(m_i c^2);
            # Calculation here is appropriate only in non-relativistic limit
            keni = self.loadDataFld('keni', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            keniint = keni - 1 - bulki
            data = keniint
        elif (text == 'bulki'):
            # Ion dimensionless bulk energy 0.5*(v_i/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
                        
            # Ion 3-velocity is zero if no particles
            v3xi[densi == 0.] = 0.
            v3yi[densi == 0.] = 0.
            v3zi[densi == 0.] = 0.
            
            # Non relativistic approximation
            bulki = 0.5*(v3xi**2 + v3yi**2 + v3zi**2)
            data = bulki
        elif (text == 'bulke'):
            # Electron dimensionless bulk energy 0.5*(v_e/c)^2;
            # Calculation here is appropriate only in non-relativistic limit
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.
            
            # Non relativistic approximation
            bulke = 0.5*(v3xe**2 + v3ye**2 + v3ze**2)
            data = bulke
        elif (text == 'vdrx'):
            # E x B drift velocity, x component
            ey = self.loadDataFld('ey', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, x component
            vdrx = (ey*bz - ez*by)/bsq

            # Cells where B = 0
            vdrx[bsq == 0.] = 0.
            data = vdrx
        elif (text == 'vdry'):
            # E x B drift velocity, y component
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            ex = self.loadDataFld('ex', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, y component
            vdry = (ez*bx - ex*bz)/bsq

            # Cells where B = 0
            vdry[bsq == 0.] = 0.
            data = vdry
        elif (text == 'vdrz'):
            # E x B drift velocity, z component
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            bsq = bx**2 + by**2 + bz**2

            # E x B / B^2, z component
            vdrz = (ex*by - ey*bx)/bsq

            # Cells where B = 0
            vdrz[bsq == 0.] = 0.
            data = vdrz
        elif (text == 'crlbx'):
            # Curl of B, x component
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            dby_dz = 0.
            dbz_dy, _ = np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])

            crlbx = dby_dz - dbz_dy
            data = crlbx
        elif (text == 'crlby'):
            # Curl of B, y component
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            _, dbz_dx= np.gradient(bz, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dbx_dz = 0.
            
            crlby = dbz_dx - dbx_dz
            data = crlby
        elif (text == 'crlbz'):
            # Curl of B, z component
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            dbx_dy, _ = np.gradient(bx, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dby_dx = np.gradient(by, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlbz = dbx_dy - dby_dx
            data = crlbz
        elif (text == 'crlex'):
            # Curl of E, x component
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            dey_dz = 0.
            dez_dy, _ = np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlex = dey_dz - dez_dy
            data = crlex
        elif (text == 'crley'):
            # Curl of E, y component
            ex = self.loadDataFld('ex', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            _, dez_dx= np.gradient(ez, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            dex_dz = 0.
            
            crley = dez_dx - dex_dz
            data = crley
        elif (text == 'crlez'):
            # Curl of E, z component
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            dex_dy, _ = np.gradient(ex, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                    1.*self.paramDict['istep']/self.paramDict['c_omp'])
            _, dey_dx= np.gradient(ey, 1.*self.paramDict['istep']/self.paramDict['c_omp'], 
                                   1.*self.paramDict['istep']/self.paramDict['c_omp'])
            
            crlez = dex_dy - dey_dx
            data = crlez
        elif (text == 'jdotE'):
            # Energization j.E, power/unit volume
            ex = self.loadDataFld('ex', self.strn1z)
            ey = self.loadDataFld('ey', self.strn1z)
            ez = self.loadDataFld('ez', self.strn1z)
            jx = self.loadDataFld('jx', self.strn1z)
            jy = self.loadDataFld('jy', self.strn1z)
            jz = self.loadDataFld('jz', self.strn1z)
            
            jdotE = ex*jx + ey*jy + ez*jz
            data = jdotE
        elif (text == 'ex_nonideal_e'):
            # Non-ideal E field, x component, electrons
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3ye[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ex_nonideal_e = ex + (v3ye*bz - v3ze*by)
            data = ex_nonideal_e
        elif (text == 'ey_nonideal_e'):
            # Non-ideal E field, y component, electrons
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3z = self.loadDataFld('v3z', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ze = ((v3z*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3zi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ze[dense == 0.] = 0.

            ey_nonideal_e = ey + (v3ze*bx - v3xe*bz)
            data = ey_nonideal_e
        elif (text == 'ez_nonideal_e'):
            # Non-ideal E field, z component, electrons
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3x = self.loadDataFld('v3x', self.strn1z)
            v3y = self.loadDataFld('v3y', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            dens = self.loadDataFld('dens', self.strn1z)
            densi = self.loadDataFld('densi', self.strn1z)
            
            # Compute electron quantities, according to density weighting
            dense = dens - densi
            v3xe = ((v3x*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3xi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            v3ye = ((v3y*(densi + (1./(self.paramDict['mi']/self.paramDict['me']))*dense) 
                     - v3yi*densi)/((1./(self.paramDict['mi']/self.paramDict['me']))*dense))
            
            # Electron 3-velocity is zero if no particles
            v3xe[dense == 0.] = 0.
            v3ye[dense == 0.] = 0.

            ez_nonideal_e = ez + (v3xe*by - v3ye*bx)
            data = ez_nonideal_e
        elif (text == 'ex_nonideal_i'):
            # Non-ideal E field, x component, ions
            ex = self.loadDataFld('ex', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ex_nonideal_i = ex + (v3yi*bz - v3zi*by)
            data = ex_nonideal_e
        elif (text == 'ey_nonideal_i'):
            # Non-ideal E field, y component, ions
            ey = self.loadDataFld('ey', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            bz = self.loadDataFld('bz', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3zi = self.loadDataFld('v3zi', self.strn1z)

            ey_nonideal_i = ey + (v3zi*bx - v3xi*bz)
            data = ey_nonideal_i
        elif (text == 'ez_nonideal_i'):
            # Non-ideal E field, z component, ions
            ez = self.loadDataFld('ez', self.strn1z)
            bx = self.loadDataFld('bx', self.strn1z)
            by = self.loadDataFld('by', self.strn1z)
            v3xi = self.loadDataFld('v3xi', self.strn1z)
            v3yi = self.loadDataFld('v3yi', self.strn1z)
            
            ez_nonideal_i = ez + (v3xi*by - v3yi*bx)
            data = ez_nonideal_i
        else:
            data = self.loadDataFld(text, self.fileRange[-1])
        
        # Choose colorbar range base on the quantity to be plotted
        lenData = len(data.flatten())
        if (text == 'dens' or text == 'densi' or text == 'dense'
            or text == 'ken' or text == 'keni'
            or text == 'keneint' or text == 'keniint'
            or text == 'bulke' or text == 'bulki'
            or text == 'pdens' or text == 'pdensi'
            or text == 'tmpxx' or text == 'tmpyy' or text == 'tmpzz'
            or text == 'tmpxxi' or text == 'tmpyyi' or text == 'tmpzzi' ):
            strt = int(lenData/100)
            zmin = 0
            zmax = np.median(np.sort(np.abs(data.flatten()))[-strt:])
        elif (text == 'bdens' or text == 'bdensi' ):
            zmin = 0
            zmax = np.max(data)
        elif (text == 'bx' or text == 'by' or text == 'bz'
              or text == 'ex' or text == 'ey' or text == 'ez'
              or text == 'ex_nonideal_e' or text == 'ey_nonideal_e' or text == 'ez_nonideal_e'
              or text == 'ex_nonideal_i' or text == 'ey_nonideal_i' or text == 'ez_nonideal_i'           
              or text == 'crlbx' or text == 'crlby' or text == 'crlbz'
              or text == 'crlex' or text == 'crley' or text == 'crlez'
              or text == 'jx' or text == 'jy' or text == 'jz'
              or text == 'v3x' or text == 'v3y' or text == 'v3z'
              or text == 'v3xi' or text == 'v3yi' or text == 'v3zi'
              or text == 'v3xe' or text == 'v3ye' or text == 'v3ze'
              or text == 'v4x' or text == 'v4y' or text == 'v4z'
              or text == 'v4xi' or text == 'v4yi' or text == 'v4zi'
              or text == 'vdrx' or text == 'vdry' or text == 'vdrz'
              or text == 'tmpxy' or text == 'tmpxz' or text == 'tmpyz'
              or text == 'tmpxyi' or text == 'tmpxzi' or text == 'tmpyzi'
              or text == 'jdotE'):
            y, x = np.histogram(data.flatten(), bins=np.min([500, len(data.flatten())]))
            ymax = np.max(y)
            indyMax = np.where(y == ymax)[0]
            if len(indyMax > 1):
                indyMax = indyMax[0]
                
            ind = np.max(np.where((y<0.05*ymax) & (x[:-1] < x[indyMax]))[0])
            zmin = x[ind]
            zmax = -1 * zmin            
        
        return zmin, zmax
        
    def getCbarRange(self):
        """
        Compute the range of the cbar (and histogram) for the type of fld file
        """
        pass
        
if __name__ == '__main__':
    app = QApplication(sys.argv)

    main = Wid()
    main.show()

    sys.exit(app.exec_())
