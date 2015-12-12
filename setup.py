'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY;
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys
from distutils.core import setup
from distutils      import sysconfig
#from distutils.extension import Extension
from numpy.distutils.core import Extension, setup

#       SETUP
#
#       12/01/2009      --      initial coding              :       jc
#       11/05/2015      --      adapted for sasmol coding   :       jc
#
#LC      1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#                                                                      *      **
'''
      setup.py is the script to install and/or update sasmol

	> sudo python setup.py install

'''

def log(logfile,line):

    print line

    logfile.write("%s\n" % (line))

    return

def compile_extensions(logfile,current_path,python):

    spath = current_path
    path = current_path+'src/extensions/'
    dpath = path+'dcdio/'
    mpath = path+'mask/'
    vpath = path+'sasview/'
    mmpath = path+'matrix_math/'

    os.chdir(dpath)
    buildst1 = python +' setup_dcdio.py build'
    result = os.popen(buildst1).readlines()
    for line in result: log(logfile,line)

    os.chdir(mpath)
    buildst2 = python +' setup_mask.py build'
    result = os.popen(buildst2).readlines()
    for line in result: log(logfile,line)

    os.chdir(vpath)
    buildst3 = python +' setup_sasview_vmd.py build'
    result = os.popen(buildst3).readlines()
    for line in result: log(logfile,line)

    os.chdir(mmpath)
    buildst4 = python +' setup_matrix_multiply.py build'
    result = os.popen(buildst4).readlines()
    for line in result: log(logfile,line)

    os.chdir(current_path)

    cpst1 = 'cp '+dpath+'build/lib*/_dcdio.so '+spath
    cpst2 = 'cp '+dpath+'dcdio.py '+spath
    cpst3 = 'cp '+mpath+'build/lib*/_mask.so '+spath
    cpst4 = 'cp '+mpath+'mask.py '+spath
    cpst5 = 'cp '+vpath+'build/lib*/_sasview_vmd.so '+spath
    cpst6 = 'cp '+vpath+'sasview_vmd.py '+spath
    cpst7 = 'cp '+mmpath+'build/lib*/matrix_math.so '+spath

    result = os.popen(cpst1).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst2).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst3).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst4).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst5).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst6).readlines()
    for line in result: log(logfile,line)
    result = os.popen(cpst7).readlines()
    for line in result: log(logfile,line)

    return

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.

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
	url='www.smallangles.net/sassie',
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

	package_dir={'src':''},

    packages=['sasmol','test_sasmol','test_sasmol/util','sasmol/extensions','sasmol/extensions/dcdio','sasmol/extensions/sasview','sasmol/extensions/mask','sasmol/extensions/matrix_math'],

	ext_modules=[
	Extension('sasmol._dcdio',['src/extensions/dcdio/dcdio.i','src/extensions/dcdio/dcdio.c'],include_dirs=[numpy_include]),
	Extension('sasmol._sasview_vmd',['src/extensions/sasview/sasview_vmd.i','src/extensions/sasview/sasview_vmd.c','src/extensions/sasview/imd.c','src/extensions/sasview/vmdsock.c'],include_dirs=[numpy_include]),
	Extension('sasmol._mask',['src/extensions/mask/mask.i','src/extensions/mask/mask.c'],include_dirs=[numpy_include]),
	Extension('sasmol.foverlap',['src/extensions/overlap/foverlap.f'],include_dirs=[numpy_include]),
	Extension('sasmol.matrix_math',['src/extensions/matrix_math/matrix_math.f'],include_dirs=[numpy_include])],
	data_files=[('src/extensions/dcdio',['src/extensions/dcdio/dcdio.i','src/extensions/dcdio/numpy.i']),('src/extensions/mask',['src/extensions/mask/mask.i','src/extensions/mask/numpy.i'])
]
	)

print '>>> copying files for installation\n\n'

current_path = os.getcwd()
spath = os.getcwd()+os.path.sep
sys.path.append('./')
python = sys.executable
logfile = open('log_sasmol_setup_install.txt','w')

print 'current_path = ',current_path
print 'current_path = ',current_path
print 'current_path = ',current_path

compile_extensions(logfile,current_path,python)
	
if(not os.path.isfile(spath+'_dcdio.so')):
	print "_dcdio.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'dcdio.py')):
	print "dcdio.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'_mask.so')):
	print "_mask.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'mask.py')):
	print "mask.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'_sasview_vmd.so')):
	print "_sasview_vmd.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'sasview_vmd.py')):
	print "sasview_vmd.py did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
if(not os.path.isfile(spath+'matrix_math.so')):
	print "matrix_math.so did not install correctly!\n\nINSTALLATION STOPPING NOW\n\n"
	sys.exit()
