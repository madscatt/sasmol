'''
    SASMOL  Copyright (C) 2009-2016 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys
from distutils.core import setup
from distutils      import sysconfig
from numpy.distutils.core import Extension, setup

#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/05/2015      --      for sasmol distribution     :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sasmol

	> sudo python setup.py install

'''

import numpy

try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(name='sasmol',
	version='2.0',
	author='Joseph E. Curtis',
	author_email='joseph.curtis@nist.gov',
	license='GPL 3',
	url='https://github.com/madscatt/sasmol',
	platforms='Linux, Mac OS X',
	description=("A library of methods to write software for molecular systems"),
	long_description=read('README.md'),	
	classifiers=["Development Status :: 4 - Beta",
		"License :: OSI Approved :: GNU Public License 3",
		"Intended Audience :: Science/Research",
		"Natural Language :: English",
		"Operating System :: Linux :: MacOS :: MacOS X",
		"Programming Language :: Python :: C :: Fortran",
		"Topic :: Scientific/Engineering :: Chemistry :: Physics"],

	package_dir={'sasmol':os.path.join('src','python')},

    packages=['sasmol','sasmol.extensions','sasmol.extensions.dcdio','sasmol.extensions.sasview','sasmol.extensions.mask','sasmol.extensions.matrix_math'],

	ext_modules=[
	Extension('sasmol._dcdio',[os.path.join('src','python','extensions','dcdio','dcdio.i'),os.path.join('src','python','extensions','dcdio','dcdio.c')],include_dirs=[numpy_include]),
	Extension('sasmol._sasview_vmd',[os.path.join('src','python','extensions','sasview','sasview_vmd.i'),os.path.join('src','python','extensions','sasview','sasview_vmd.c'),os.path.join('src','python','extensions','sasview','imd.c'),os.path.join('src','python','extensions','sasview','vmdsock.c')],include_dirs=[numpy_include]),
	Extension('sasmol._mask',[os.path.join('src','python','extensions','mask','mask.i'),os.path.join('src','python','extensions','mask','mask.c')],include_dirs=[numpy_include]),
	Extension('sasmol.foverlap',[os.path.join('src','python','extensions','overlap','foverlap.f')],include_dirs=[numpy_include]),
	Extension('sasmol.matrix_math',[os.path.join('src','python','extensions','matrix_math','matrix_math.f')],include_dirs=[numpy_include])],
	data_files=[(os.path.join('src','python','extensions','dcdio'),[os.path.join('src','python','extensions','dcdio','dcdio.i'),os.path.join('src','python','extensions','dcdio','numpy.i')]),(os.path.join('src','python','extensions','mask'),[os.path.join('src','python','extensions','mask','mask.i'),os.path.join('src','python','extensions','mask','numpy.i')])
]
	)


