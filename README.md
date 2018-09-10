Widpy
=====

Widpy is a (python 3.6, PyQt5, PyQtgraph) interactive tool for visualization of
data from the particle-in-cell code TRISTAN-MP, in the same spirit as the IDL
widget `wid.pro` (authors: Anatoly Spitkovsky, Lorenzo Sironi, Uri Keshet) and python
package  `Iseult` (author: Patrick Crumley).  Here is a short Widpy demonstration:

![widpy_demo_480p](https://user-images.githubusercontent.com/38045958/45279979-b3da5f00-b4a0-11e8-9ecf-3708ee691d06.gif)

While intended for general-purpose visualization of TRISTAN-MP data, Widpy
includes features that are specific to simulations of magnetic reconnection:

* Option to display a contour delineating the reconnection region, according to user selected threshold value
* Slider controlling the reconnection threshold value, `dthresh`; the reconnection region is selected based on a mixing criterion between right-tagged and left-tagged particles.  Cells are counted as ones where reconnection has occurred if the mixing (between left- and right-tagged particles) criterion

```
dthresh < (right-tagged particle density)/(total density) < 1-dthresh
```

is satisfied.

Other features:

* ROI (region of interest) mode for plotting spectra of particles
* Several matplotlib colorbars (viridis, RdBu, plasma, cubehelix)

Widpy is tested with Linux Red Hat, but no guarantees it will work with other operating systems.  Development is ongoing, but feel free to submit an [issue](https://github.com/mrowan137/widpy/issues).

Copyright 2018 Michael E. Rowan, Cambridge MA


Prerequisites
-------------
To use Widpy:

* Python 3.6.x
* NumPy
* Matplotlib
* Pyqt
* Pyqtgraph
* H5py

and to run tests:

* Nose

It is recommended to install these using a package manager, e.g. Conda:

```
conda create -n py36 python=3.6
```

and then install each package listed above:

```
conda install numpy
conda install matplotlib
...
```

See [requirements.txt](requirements.txt) for details.


Installation
------------

If the necessary packages are installed, simply clone the repository and navigate to a directory containing tristan `/output` data (or alternatively, the `/output` directory; either should work).  Widpy can be run with:

```
python /path/to/widpy/widpy/widpy.py
```

Widpy is intended to be run from the command line; to make it globally executable, navigate to `/path/to/widpy/widpy` and run

```
chmod +x widpy.py
```
    
and add to your `PATH` variable, e.g.

```
PATH=/path/to/widpy/widpy/:$PATH
export PATH
```

To execute Widpy, navigate to a directory containing TRISTAN data and run:

```
widpy.py
```

Running the tests
-----------------

To test, navigate to `/path/to/widpy/tests` and run

```
python -m nose
```

Usage
-----

Summary of the controls (graphic below):

* Navigation bar - type a directory containing TRISTAN data and press enter to load data from that location
* Dropdown menus - select a quantity to plot; top, middle, and bottom dropdown menus correspond to top, middle, and bottom displayed images
* Output file slider - click and drag slider to select file number
* Pushbuttons  - linked to the slider; use to step through TRISTAN data according to the step size (displayed between pushbuttons; the step size is editable)
* FPS slider  - frames per second for video mode
* Rec. slider - set the threshold value `dthresh` for determination of reconnected cells
* Prtl downsample slider - downsampling factor for particle files; this can be used to speed up display of particle spectra in the ROI (region of interest)
* Parameter list - displays parameters read from the directory input file
* ROI - 'Region of Interest'; activate to display a (draggable, resizable, rotatable) box on the image, and a second plot window to display one of the following:
  * 1D average (along the zeroth axis of the box) of data in the ROI
  * Log-log histogram of particle number vs. gamma-1 of particles in the ROI
* Menu - dropdown menu with several options:
  * ROI mode - toggles display quantity when ROI is actiated; if checked, ROI shows log-log histogram of particle number vs. gamma-1, if not checked ROI shows 1D average of data within the ROI
  * Rec region - toggles display of isocontour delineating cells where reconnection has occurred, according to the mixing criterion referenced above; threshold value `dthresh` is selected with dthresh slider on the left
  * Matplotlib cbar - select a matplotlib colorbar (viridis, RdBu, plasma, cubehelix)

The ImageView object (consisting of the image object, and the colorbar/histogram) is from the PyQtGraph library and contains several builtin controls:
* scroll up/down [on image] - zoom out/in
* right click [on image] - show display and export options
* spacebar [on image] - 'play video'; rate is controlled by FPS slider, and step size controlled by value between pushbuttons
* left/right [on image] - same function as the pushbuttons; loop through files, with step size as displayed between pushbuttons
* right click [on histogram] - access axes control for histogram
* scroll up/down [on histogram] - increase/decrease the histogram range; drag the lines extending from the colorbar to adjust the mapping
* right click [on colorbar] - select a builtin colorbar, as well as other colorbar options

Read more about the ImageView object in the [PyQtgraph documentation](http://pyqtgraph.org/documentation/)

The basic controls referenced above are illustrated in this graphic:

![widpy_tutorial 001](https://user-images.githubusercontent.com/38045958/45284161-cd35d800-b4ad-11e8-82ca-7ac3a69906cc.jpeg)


Built With
----------

* [Anaconda](https://docs.anaconda.com/)
* [H5py](http://docs.h5py.org/en/latest/) 
* [Matplotlib](https://matplotlib.org/contents.html)
* [NumPy](https://docs.scipy.org/doc/numpy/)
* [PyQt5](http://pyqt.sourceforge.net/Docs/PyQt5/)
* [PyQtgraph](http://pyqtgraph.org/documentation/)


License
-------

This project is licensed under the GPLv3 License - see the [LICENSE.md](LICENSE.md) file for details.


Acknowledgments
---------------

* Anatoly Spitkovsky, Lorenzo Sironi, Uri Keshet, authors of IDL widget wid.pro
* Patrick Crumley, author of python package [Iseult](https://github.com/pcrumley/Iseult)
* Luke Campagnola, author of [PyQtgraph](http://pyqtgraph.org/documentation/)

