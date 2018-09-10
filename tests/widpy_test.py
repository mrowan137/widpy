"""
widpy_test.py - Test module for visualization widget, using sample TRISTAN data
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
from PyQt5.QtTest import QTest

# Qt based plotting widgets
import pyqtgraph as pg
#from pyqtgraph import ptime

# For HDF5 files
import h5py

# Main widget
from widpy.widpy import Wid

# For unit test
import unittest
from unittest import TestCase

app = QApplication(sys.argv)
class WidpyTest(TestCase):
    """
    Test of the main widget; initialize with sample data
    """
    def testInstance(self):
        """
        Instantiate widget
        """
        wid = Wid()
        #assert(wid.fileRoot ==)
                
    def testLineEdit(self):
        """
        Check default fileroot
        """
        wid = Wid()
        defaultRoot = (os.path.dirname(os.getcwd()) + '/widpy/../tests/test_data/')
        print(defaultRoot)
        print(str(wid.rootLine.text()))
        assert(str(wid.rootLine.text()) == defaultRoot)
        
    def testSliderInit(self):
        """
        Check init values for sliders
        """
        wid = Wid()
        assert(wid.sl.value() == 1)
        assert(wid.slRec.value() == 30)
        assert(wid.slPrtlDown.value() == 10)
        
    def testComboBoxInit(self):
        """
        Test combo box init strings
        """
        wid = Wid()
        assert(str(wid.comboBox1.currentText()) == 'dens')
        assert(str(wid.comboBox2.currentText()) == 'v3y')
        assert(str(wid.comboBox3.currentText()) == 'jz')
        
    def testPushButtons(self):
        """
        Test push-button/slider interface; click left, right and check slider value
        """
        wid = Wid()
        QTest.mouseClick(wid.TB1, Qt.LeftButton)
        assert(wid.sl.value() == 5)
        
        QTest.mouseClick(wid.TB2, Qt.LeftButton)
        assert(wid.sl.value() == 1)
        
    def testImageView(self):
        """
        Test imageView/slider interface; click right, left, right for imv1, imv2,
        imv3, and check the slider value at each step
        """
        wid = Wid()
        QTest.keyClick(wid.imv1, Qt.Key_Right)
        assert(wid.sl.value() == 2)
        QTest.keyClick(wid.imv2, Qt.Key_Left)
        assert(wid.sl.value() == 1)
        QTest.keyClick(wid.imv3, Qt.Key_Right)
        assert(wid.sl.value() == 2)
        
    
if __name__ == '__main__':
    unittest.main()
