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

#__init__.py     test_sascalc/       test_sasop/
#data/           test_sasio/     test_sasproperties/
#manual_tests/       test_sasmath/       test_sassubset/
#script/         test_sasmol/        util/

    packages=['sasmol','sasmol.test_sasmol','sasmol.test_sasmol.util','sasmol.test_sasmol.manual_tests','sasmol.test_sasmol.data','sasmol.test_sasmol.data.pdb_common','sasmol.test_sasmol.data.dcd_common','sasmol.test_sasmol.data.sasmol','sasmol.test_sasmol.data.sasmol.sascalc','sasmol.test_sasmol.data.sasmol.sasio','sasmol.test_sasmol.data.sasmol.sasmath','sasmol.test_sasmol.data.sasmol.sasmol','sasmol.test_sasmol.data.sasmol.sasop','sasmol.test_sasmol.data.sasmol.sasproperties','sasmol.test_sasmol.test_sascalc','sasmol.test_sasmol.test_sasio','sasmol.test_sasmol.test_sasmath','sasmol.test_sasmol.test_sasmol','sasmol.test_sasmol.test_sasop','sasmol.test_sasmol.test_sasproperties','sasmol.test_sasmol.test_sassubset','sasmol.extensions','sasmol.extensions.dcdio','sasmol.extensions.sasview','sasmol.extensions.mask','sasmol.extensions.matrix_math'],
	
    ext_modules=[
	Extension('sasmol._dcdio',[os.path.join('src','python','extensions','dcdio','dcdio.i'),os.path.join('src','python','extensions','dcdio','dcdio.c')],include_dirs=[numpy_include]),
	Extension('sasmol._sasview_vmd',[os.path.join('src','python','extensions','sasview','sasview_vmd.i'),os.path.join('src','python','extensions','sasview','sasview_vmd.c'),os.path.join('src','python','extensions','sasview','imd.c'),os.path.join('src','python','extensions','sasview','vmdsock.c')],include_dirs=[numpy_include]),
	Extension('sasmol._mask',[os.path.join('src','python','extensions','mask','mask.i'),os.path.join('src','python','extensions','mask','mask.c')],include_dirs=[numpy_include]),
	Extension('sasmol.foverlap',[os.path.join('src','python','extensions','overlap','foverlap.f')],include_dirs=[numpy_include]),
	Extension('sasmol.matrix_math',[os.path.join('src','python','extensions','matrix_math','matrix_math.f')],include_dirs=[numpy_include])],

    data_files = [ 
        ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag.pdb')]),
        ( os.path.join('sasmol','test_sasmol','manual_tests') , [os.path.join('src','python','test_sasmol','manual_tests','hiv1_gag_200_frames.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','1ATM.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','2AAD.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','rna-1to10.dcd')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t5.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t6.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','dcd_common') , [os.path.join('src','python','test_sasmol','data','dcd_common','t7.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1ATM-1to2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1ATM.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1KP8.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','2AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','rna.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1PSI.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','3AAD-2chain.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','dimcd_fixed_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','1CRN.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','2AAD-1to3.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','3AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','pdb_common') , [os.path.join('src','python','test_sasmol','data','pdb_common','rna-1to10.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','sascalc') , [os.path.join('src','python','test_sasmol','data','sasmol','sascalc','1ATN.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sascalc') , [os.path.join('src','python','test_sasmol','data','sasmol','sascalc','1CRN-rot.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sascalc') , [os.path.join('src','python','test_sasmol','data','sasmol','sascalc','1CRN-rot-shift.pdb')]),
    
    
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Catoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','ConflictAtoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Hatoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','MisAtoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Natoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Oatoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Otheratoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Patoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','Satoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','amino_acid_sld.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','charmm27_atoms.txt')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasproperties') , [os.path.join('src','python','test_sasmol','data','sasmol','sasproperties','standard_atomic_weigh.txt')]),
   
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasop') , [os.path.join('src','python','test_sasmol','data','sasmol','sasop','1CRN-rot-shift.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasop') , [os.path.join('src','python','test_sasmol','data','sasmol','sasop','1CRN-rot-sub.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasop') , [os.path.join('src','python','test_sasmol','data','sasmol','sasop','1CRN-rot.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasop') , [os.path.join('src','python','test_sasmol','data','sasmol','sasop','1CRN-sub.pdb')]),
    
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasmath') , [os.path.join('src','python','test_sasmol','data','sasmol','sasmath','1CRN-rot-shift.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasmath') , [os.path.join('src','python','test_sasmol','data','sasmol','sasmath','1CRN-rot.pdb')]),
        
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasmol') , [os.path.join('src','python','test_sasmol','data','sasmol','sasmol','1CRN-3frames.pdb')]),

        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','1AA-NoEND.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','1ATM-1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','1ATM-2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-END_wrong_number_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-MODEL.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-MODEL_missing_END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-MODEL_mix_END_noterminating.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-MODEL_wrong_number_atoms.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3-MODEL_wrongnumber_mix_END.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3_MODEL.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-1to3_MODELwrong.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-2.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','2AAD-3.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','3AAD-2chain.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','3AAD.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','nef_nohis.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','nef_nohis_1.pdb')]),
        ( os.path.join('sasmol','test_sasmol','data','sasmol','sasio') , [os.path.join('src','python','test_sasmol','data','sasmol','sasio','new_package_rna.pdb')]) #,

                   ]
	)

# (jc@air)sasmol/sascalc% ls
#1ATN.pdb        1CRN-rot-shift.pdb  1CRN-rot.pdb

