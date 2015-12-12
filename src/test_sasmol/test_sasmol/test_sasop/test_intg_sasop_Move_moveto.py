'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
	 Core-Testing: Copyright (C) 2011 Hailiang Zhang, Ph.D.

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

from sassie.core_testing.util import env,util



from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS
from sassie.sasmol import sasmol, sasop, sascalc

import numpy


import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'

class Test_intg_sasop_Move_moveto(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)

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
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        with self.assertRaises(Exception):
          self.o.moveto(0,value)

    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        result_com = self.o.com()
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        expected_com = value
        self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)


    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        result_com = self.o.com()
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        expected_com = value
        #self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)

    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        result_com = self.o.com()
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        expected_com = value
        #self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)

    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        result_com = self.o.com()
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        expected_com = value
        #self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        result_com = self.o.com()
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        expected_com = value
        #self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

