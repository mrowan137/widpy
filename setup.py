from setuptools import setup
import os

setup(name='Widpy',
      version='0.1',
      description='Widpy: TRISTAN data visualization',
      url='http://github.com/mrowan137/widpy',
      author='Michael Rowan',
      author_email='mrowan137@gmail.com',
      license='GPLv3',
      packages=['widpy'],
      install_requires=['numpy',###
			'matplotlib',
			'pyqt',###
                        'pyqtgraph',
                        'h5py'],
      test_suite='nose.collector',
      tests_require=['nose'],
)
