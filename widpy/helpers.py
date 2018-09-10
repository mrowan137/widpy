"""
helpers.py - Classes/functions for use in widpy.py
"""
# Standard imports
import numpy as np

# PyQt5 Widgets
from PyQt5.QtCore import Qt

# Qt based plotting widgets
import pyqtgraph as pg
from pyqtgraph import ptime


class ImageView(pg.ImageView):
    """
    Customization of pyqtgraph.ImageView, from the pyqtgraph library by Luke Campagnola
    """
    def __init__(self, roiPos=(0,0), roiSize=(10,10), parent=None):
        super(ImageView, self).__init__()
        # ROI, Rec initialization
        self.roiMode = False
        self.recMode = False
        self.roi.setPos(roiPos)
        self.roi.setSize(roiSize)
        self.roiCurves = []
        self.scatterx = None
        self.scattery = None
        self.scatterz = None
        
    def updateImage(self, autoHistogramRange=False): #autoRange=False
        """
        Redraw image on screen, and if ROI mode selected, the scatter plot
        """
        # Redraw scatter on screen
        if self.image is None:
            return
            
        image = self.image
        
        if self.roiMode:
            self.plt.setData(self.scattery, self.scatterx)
            
        if self.recMode:
            self.iso.setData(self.parent.recReg)
        
        if autoHistogramRange:
            self.ui.histogram.setHistogramRange(self.levelMin, self.levelMax)
        
        # Transpose image into order expected by ImageItem
        if self.imageItem.axisOrder == 'col-major':
            axorder = ['t', 'x', 'y', 'c']
        else:
            axorder = ['t', 'y', 'x', 'c']
        axorder = [self.axes[ax] for ax in axorder if self.axes[ax] is not None]
        image = image.transpose(axorder)
                        
        self.imageItem.updateImage(image)
         
    def buildMenu(self):
        """
        Create menu

        Selections:
            *ROI mode - If checked, plot the particle spectra in the ROI.
                        If not checked, plot 1D average of the fluid data in ROI.
            *Rec mode - If selected, plot isocontour of (density of right-tagged
                        particles)/(total density); threshhold value is controlled
                        with `dthresh` threshold, set with slRec slider
            *matplotlib cbar - select a colortable from matplotlib 
        """
        self.menu = pg.QtGui.QMenu()
        
        # ROI mode
        self.roiAction = pg.QtGui.QAction("ROI mode", self.menu)
        self.roiAction.setCheckable(True)    
        self.roiAction.triggered.connect(self.roiToggled)
        self.menu.addAction(self.roiAction)
        
        # Rec mode
        self.recAction = pg.QtGui.QAction("Rec region", self.menu)
        self.recAction.setCheckable(True)
        self.recAction.triggered.connect(self.recToggled)
        # Include Rec option only if necessary arrays are present in output
        try:
            pdensr = self.parent.loadDataFld('pdens', self.parent.strn1z)
            bdens = self.parent.loadDataFld('bdens', self.parent.strn1z)
            recArrsPresent = True
        except KeyError:
            print('Arrays needed for identification of reconnection region not found!')
            recArrsPresent = False
            
        if (recArrsPresent):
            self.menu.addAction(self.recAction)
        else:
            pass
        
        # Matplotlib colortables
        self.cbarMenu = self.menu.addMenu('matplotlib cbar')
        
        self.viridisAction = pg.QtGui.QAction('viridis', self.cbarMenu)
        self.viridisAction.triggered.connect(self.viridisCbar)
        self.cbarMenu.addAction(self.viridisAction)
        
        self.RdBuAction = pg.QtGui.QAction('RdBu', self.cbarMenu)
        self.RdBuAction.triggered.connect(self.RdBuCbar)
        self.cbarMenu.addAction(self.RdBuAction)

        self.plasmaAction = pg.QtGui.QAction('plasma', self.cbarMenu)
        self.plasmaAction.triggered.connect(self.plasmaCbar)
        self.cbarMenu.addAction(self.plasmaAction)

        self.cubehelixAction = pg.QtGui.QAction('cubehelix', self.cbarMenu)
        self.cubehelixAction.triggered.connect(self.cubehelixCbar)
        self.cbarMenu.addAction(self.cubehelixAction)
                
        
    def roiChanged(self):
        """
        Update display of data selected by ROI
        """ 
        image = self.image

        if self.image is None:
            return
                
        # Extract image data from ROI
        axes = (self.axes['x'], self.axes['y'])
        data, coords = self.roi.getArrayRegion(image.view(np.ndarray), 
                                               self.imageItem, 
                                               axes, returnMappedCoords=True)
        
        # Extract slice from ROI
        wh = self.roi.getArraySlice(image.view(np.ndarray), self.imageItem, axes)[0]
        arr = np.zeros(np.shape(image.view(np.ndarray)))
        arr[wh] = 1.
        
        if (self.roiMode):
            # Find where scatter data overlap with the ROI
            # Particles selected if ROI overlaps with a cell containing a particle
            if self.scatterx is None:
                return
            
            my, mx = arr.shape
            whIn = np.where((((self.scatterx.astype(int) - 0) >= 0) &
                             ((self.scatterx.astype(int) - 0) < mx) &
                             ((self.scattery.astype(int) - 0) >= 0) &
                             ((self.scattery.astype(int) - 0) < my)))[0]
            xInd = (self.scatterx.astype(int) - 0)
            yInd = (self.scattery.astype(int) - 0)

            wh = np.where(arr[yInd[whIn], xInd[whIn]] == 1)[0]
            whp = (whIn[wh])
            
            if len(whp != 0):
                data_prtl = self.scatterz[whp]
            else:
                data_prtl = np.zeros(1)
                
            if data is None or data_prtl is None:
                return

            # Convert extracted data into 1D plot data
            # Bin the gammas
            log_gamma_min = np.min(np.log10(self.scatterz-1.))
            log_gamma_max = np.max(np.log10(self.scatterz-1.))
            
            nbins = 20
            
            logbins=np.logspace(log_gamma_min, log_gamma_max, nbins+1)
            logbins_center = 10**(0.5*(
                    np.log10(logbins[1::]) + np.log10(logbins[0:-1])))
            h = (np.histogram(data_prtl-1., logbins, normed=True))[0]
            
            # If all -inf slice, fill with zeros (avoid error warning)
            if np.isnan(h).all():
                h[:] = 1.
            else:
                h[h == -np.inf] = np.nan
            
            xvals = np.log10(logbins_center)
            data = np.log10(h)
        else:
            # Average across y-axis of ROI
            data = data.mean(axis=axes[1])
            coords = coords[:,:,0] - coords[:,0:1,0]
            xvals = (coords**2).sum(axis=0) ** 0.5
                            
        plots = [(xvals, data, 'w')]

        # Update plot
        while len(plots) < len(self.roiCurves):
            c = self.roiCurves.pop()
            c.scene().removeItem(c)
        while len(plots) > len(self.roiCurves):
            self.roiCurves.append(self.ui.roiPlot.plot())

        x, y, p = plots[0]
        self.roiCurves[0].setData(x, y, pen=p)

    def roiToggled(self, b):
        """
        Set the ROI mode and update image
        """
        self.roiMode = b
        if b:
            # Initialize the scatterplot
            xe = (self.parent.loadDataPrtl('xe', self.parent.strn1z)
                  /self.parent.paramDict['istep'] + self.parent.xshift)
            ye = (self.parent.loadDataPrtl('ye', self.parent.strn1z)
                  /self.parent.paramDict['istep'])
            gammae = self.parent.loadDataPrtl('gammae', self.parent.strn1z)
            
            self.setScatter((xe,ye,gammae))
            self.plt = pg.ScatterPlotItem(ye, xe, 
                                       pen=pg.mkPen(None), brush=pg.mkBrush(66, 255, 51, 100))
            self.vb = self.imageItem.getViewBox()
            self.vb.addItem(self.plt)
        if not b:
            self.vb.removeItem(self.plt)
            
        self.roiChanged()
        self.updateImage()
        
    def recToggled(self, b):
        """
        Set the rec mode and update image
        """
        self.recMode = b
        if b:
            # Enable the dthresh bar
            self.parent.slRec.setEnabled(True)
             
            # Make the isocurve plot item
            self.iso = pg.IsocurveItem(level=0.8, pen='w')
            self.vb = self.imageItem.getViewBox()
            self.vb.addItem(self.iso)
            self.iso.setZValue(1.)
            
            # Initialize the rec region
            pdensrFull = self.parent.loadDataFld('pdens', self.parent.strn1z)
            pdensr = self.parent.subsetData(pdensrFull)
            densFull = self.parent.loadDataFld('dens', self.parent.strn1z)
            dens = self.parent.subsetData(densFull)
            bdensFull = self.parent.loadDataFld('bdens', self.parent.strn1z)
            bdens = self.parent.subsetData(bdensFull) 
            
            pdenstot = dens - bdens
            
            rat = self.parent.nandiv(pdensr, pdenstot)
            whrec=np.where((rat > self.parent.dthresh) & 
                                (rat < (1-self.parent.dthresh)) &
                                (bdens == 0))
            
            my, mx = np.shape(bdens)
            self.parent.recReg = np.zeros((my, mx))
            self.parent.recReg[whrec]=1.
        if not b:
            self.parent.slRec.setEnabled(False)
            self.vb.removeItem(self.iso)
            
        self.updateImage()
                    
    def setScatter(self, data):
        """
        Set data (x, y, z) for scatter
        """
        self.scatterx = data[0]
        self.scattery = data[1]
        self.scatterz = data[2]
        
    def play(self, key):
        """
        Deactivate the default pg.ImageView play
        """
        pass
    
    def keyPressEvent(self, ev):
        """
        Control of frame display using keys
            *space - play movie at the rate FPS 
            *left - frame --> frame - 1
            *right - frame --> frame + 1
        """
        if ev.key() == Qt.Key_Space:
            if self.playRate == 0:
                self.parent.slFPS.setEnabled(False)
                fps = self.parent.fps
                self.onVideo(fps)
                ev.accept()
            else:
                self.parent.slFPS.setEnabled(True)
                self.onVideo(0)
                ev.accept()
        elif ev.key() == Qt.Key_Left:
            self.parent.onActivateTB1()
            ev.accept()
        elif ev.key() == Qt.Key_Right:
            self.parent.onActivateTB2()
            ev.accept()
                
    def onVideo(self, rate):
        """
        Loop through frames at FPS = `rate`
        """
        self.playRate = rate
        if rate == 0:
            self.playTimer.stop()
            return
            
        self.lastPlayTime = ptime.time()
        if not self.playTimer.isActive():
            # 1 ms timer; emit PySide.QtCore.QTimer.timeout() at 1ms interval
            self.playTimer.start(1)

    def jumpFrames(self, _):
        """
        Override the default ImageView jump frame, link to frame jump functions
        """
        if (self.parent.fps > 0):
            self.parent.onActivateTB2()
        elif (self.parent.fps <0):
            self.parent.onActivateTB1()
            
    def viridisCbar(self):
        """
        Use Matplotlib colortable `viridis` for ImageView colorbar
        """
        self.parent.cbarName = 'viridis'
        
        # Initialize a matplotlib colortable
        self.hist = self.getHistogramWidget()
        self.parent.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.parent.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist.gradient.ticks:
            tick.hide()
            
    def RdBuCbar(self):
        """
        Use Matplotlib colortable `RdBu` for ImageView colorbar
        """
        self.parent.cbarName = 'RdBu'
        
        # Initialize a matplotlib colortable
        self.hist = self.getHistogramWidget()
        self.parent.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.parent.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist.gradient.ticks:
            tick.hide()
            
    def plasmaCbar(self):
        """
        Use Matplotlib colortable `plasma` for ImageView colorbar
        """
        self.parent.cbarName = 'plasma'
        
        # Initialize a matplotlib colortable
        self.hist = self.getHistogramWidget()
        self.parent.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.parent.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist.gradient.ticks:
            tick.hide()
            
    def cubehelixCbar(self):
        """
        Use Matplotlib colortable `cubehelix` for ImageView colorbar
        """
        self.parent.cbarName = 'cubehelix'
        
        # Initialize a matplotlib colortable
        self.hist = self.getHistogramWidget()
        self.parent.setCbar()
        pos, color = np.arange(0, 1, 1/256), self.parent.colorTable
        
        # Hide the ticks on colorbar (only need for matplotlib cmap)
        self.hist.gradient.setColorMap(pg.ColorMap(pos, color))
        for tick in self.hist.gradient.ticks:
            tick.hide()
            
