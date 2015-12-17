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

from sasmol.test_sasmol.util import env, util

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol

import numpy

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sascalc_Prop_calcminmax(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)


    def test_null(self):
        with self.assertRaises(Exception):
            self.o.calcminmax_frame(0)


    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        result_minmax  = self.o.calcminmax_frame(0)
        expected_minmax = [ [73.944, 41.799, 41.652], [73.944, 41.799, 41.652]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        result_minmax  = self.o.calcminmax_frame(0)
        print result_minmax
        expected_minmax = [ [70.721,  41.799,  39.354],[ 79.712,  46.273,  43.910]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_two_aa_3frames_pdb(self):
        self.o.read_pdb(DataPath+'2AAD-1to3.pdb')
        result_minmax  = self.o.calcminmax_frame(2)
        print result_minmax
        expected_minmax = [[70.721, -46.273, 39.354], [ 79.712, -41.799, 43.910]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calcminmax_frame(0)
        print result_minmax
        expected_minmax = [[-88.148, -86.246, -81.494],[ 84.491, 77.158, 84.429]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        result_minmax  = self.o.calcminmax_frame(0)
        print result_minmax
        expected_minmax = [[-3.097, -0.516, -7.422],[24.284, 20.937, 19.580]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        result_minmax  = self.o.calcminmax_frame(0)
        print result_minmax
        expected_minmax = [[8.043, -73.261, -48.819], [156.999, 75.260, 101.562]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

