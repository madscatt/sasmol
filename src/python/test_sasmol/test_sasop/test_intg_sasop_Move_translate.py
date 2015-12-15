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

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_intg_sasop_Move_translate(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)
        self.tol = 3

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
          self.o.translate(0,value)


    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        print result_coor
        expected_coor = numpy.array([[[74.944, 44.799, 47.652]]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,self.tol)



    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor = self.o.coor()
        numpy.set_printoptions(precision=3)
        print result_coor
        expected_coor = numpy.array([[[ 74.944, 44.799, 47.652],\
                                      [ 75.229, 45.563, 46.456],\
                                      [ 76.667, 46.093, 46.463],\
                                      [ 77.264, 46.279, 45.401],\
                                      [ 74.210, 46.734, 46.336],\
                                      [ 72.856, 46.168, 45.926],\
                                      [ 74.677, 47.782, 45.354],\
                                      [ 71.721, 47.177, 45.946],\
                                      [ 77.231, 46.330, 47.647],\
                                      [ 78.592, 46.852, 47.730],\
                                      [ 79.617, 45.820, 48.184],\
                                      [ 80.712, 46.169, 48.656],\
                                      [ 78.671, 48.097, 48.648],\
                                      [ 78.054, 47.816, 49.910],\
                                      [ 77.970, 49.273, 48.000]]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,self.tol)


    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor_sum = sum(sum(sum(self.o.coor())))
        expected_coor_sum = 165825.827
        self.assertAlmostEqual(expected_coor_sum, result_coor_sum,self.tol)


    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor_sum = sum(sum(sum(self.o.coor())))
        expected_coor_sum = 11779.587
        self.assertAlmostEqual(expected_coor_sum, result_coor_sum,self.tol)


    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.translate(0,value)
        result_coor_sum = sum(sum(sum(self.o.coor())))
        expected_coor_sum = 6840020.260
        self.assertAlmostEqual(expected_coor_sum, result_coor_sum,self.tol)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

