'''
    SASMOL: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

from sasmol.test_sasmol.util import env,util

from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop
import sasmol.sascalc as sascalc

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasop')+os.path.sep

class Test_intg_sasop_Move_align(MockerTestCase): 


    def setUp(self):
        self.o1=sasmol.SasMol(0)
        self.o2=sasmol.SasMol(0)
        self.o1Sub=sasmol.SasMol(0)
        self.o2Sub=sasmol.SasMol(0)


    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if isinstance(a[i],(int,float,numpy.generic)):
                 if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
                 self.assertAlmostEqual(a[i],b[i],places)
              else:
                 self.assert_list_almost_equal(a[i],b[i],places)


    def test_null(self):
        with self.assertRaises(Exception):
            self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)


    def test_1ATM_against_1ATMRot_subset_full_pdb(self):
        self.o1.read_pdb(DataPath+'1ATM.pdb')
        self.o2.read_pdb(DataPath+'1ATM.pdb')
        self.o1Sub.read_pdb(DataPath+'1ATM.pdb')
        self.o2Sub.read_pdb(DataPath+'1ATM.pdb')

        frame = 0
        comO1Sub = self.o1Sub.calccom(frame)
        comO2Sub = self.o2Sub.calccom(frame)
        self.o1Sub.center(frame)
        self.o2Sub.center(frame)
        coorO1Sub = self.o1Sub.coor()[frame]
        coorO2Sub = self.o2Sub.coor()[frame]

        expected_coor = self.o1.coor()[frame]
        self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)

        result_coor = self.o1.coor()[frame]
        print 'expected:\n', expected_coor
        print 'result:\n', result_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_2AAD_against_2AADRot_subset_full_pdb(self):
        self.o1.read_pdb(DataPath+'2AAD.pdb')
        self.o2.read_pdb(DataPath+'2AAD.pdb')
        self.o1Sub.read_pdb(DataPath+'2AAD.pdb')
        self.o2Sub.read_pdb(DataPath+'2AAD.pdb')

        frame = 0
        comO1Sub = self.o1Sub.calccom(frame)
        comO2Sub = self.o2Sub.calccom(frame)
        self.o1Sub.center(frame)
        self.o2Sub.center(frame)
        coorO1Sub = self.o1Sub.coor()[frame]
        coorO2Sub = self.o2Sub.coor()[frame]

        expected_coor = self.o1.coor()[frame]
        self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)

        result_coor = self.o1.coor()[frame]
        print 'expected:\n', expected_coor
        print 'result:\n', result_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_1CRN_against_1CRNRot_subset_full_pdb(self):
        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(moduleDataPath+'1CRN-rot.pdb')
        self.o1Sub.read_pdb(DataPath+'1CRN.pdb')
        self.o2Sub.read_pdb(moduleDataPath+'1CRN-rot.pdb')

        frame = 0
        comO1Sub = self.o1Sub.calccom(frame)
        comO2Sub = self.o2Sub.calccom(frame)
        self.o1Sub.center(frame)
        self.o2Sub.center(frame)
        coorO1Sub = self.o1Sub.coor()[frame]
        coorO2Sub = self.o2Sub.coor()[frame]

        expected_coor = self.o1.coor()[frame]
        self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)

        result_coor = self.o1.coor()[frame]
        print 'expected:\n', expected_coor
        print 'result:\n', result_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_1CRN_against_1CRNRotShift_subset_full_pdb(self):
        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(moduleDataPath+'1CRN-rot-shift.pdb')
        self.o1Sub.read_pdb(DataPath+'1CRN.pdb')
        self.o2Sub.read_pdb(moduleDataPath+'1CRN-rot-shift.pdb')

        frame = 0
        comO1Sub = self.o1Sub.calccom(frame)
        comO2Sub = self.o2Sub.calccom(frame)
        self.o1Sub.center(frame)
        self.o2Sub.center(frame)
        coorO1Sub = self.o1Sub.coor()[frame]
        coorO2Sub = self.o2Sub.coor()[frame]

        expected_coor = self.o1.coor()[frame]
        self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)

        result_coor = self.o1.coor()[frame]
        print 'expected:\n', expected_coor
        print 'result:\n', result_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_1CRN_against_1CRNRotShift_subset_aa1to3_pdb(self):
        self.o1.read_pdb(DataPath+'1CRN.pdb')
        self.o2.read_pdb(moduleDataPath+'1CRN-rot-shift.pdb')
        self.o1Sub.read_pdb(moduleDataPath+'1CRN-sub.pdb')
        self.o2Sub.read_pdb(moduleDataPath+'1CRN-rot-sub.pdb')

        frame = 0
        comO1Sub = self.o1Sub.calccom(frame)
        comO2Sub = self.o2Sub.calccom(frame)
        self.o1Sub.center(frame)
        self.o2Sub.center(frame)
        coorO1Sub = self.o1Sub.coor()[frame]
        coorO2Sub = self.o2Sub.coor()[frame]

        expected_coor = self.o1.coor()[frame]
        self.o1.align(frame, coorO2Sub, comO2Sub, coorO1Sub, comO1Sub)

        result_coor = self.o1.coor()[frame]
        print 'expected:\n', expected_coor
        print 'result:\n', result_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)



    def tearDown(self):
        pass


if __name__ == '__main__': 
   main() 

