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
dcdDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep

class Test_sascalc_Prop_calc_minmax_all_steps(MockerTestCase): 

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


    def test_one_atom(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'1ATM.dcd')
        expected_minmax = [ [73.944, 38.799, 41.652], [76.944,  41.799, 41.652]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_two_aa(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'2AAD.dcd')
	print result_minmax
        expected_minmax = [[-79.712, -46.273,  39.354], [79.712,  46.273,  43.910]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])


    def test_rna_1to10(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-1to10.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-43.801, -44.888, -42.605], [ 41.234,  39.706,  41.903]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
    def test_rna_0point8g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-0.8g.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-45.714, -45.643, -42.868], [ 41.65 ,  41.087,  45.362]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
    def test_rna_1point0g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-1.0g.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-46.369, -45.643, -42.868], [ 41.65 ,  41.087,  45.362]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
    def test_rna_2point0g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-2.0g.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-48.253, -45.643, -42.868], [ 41.65 ,  41.087,  45.362]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
    def test_rna_3point2g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-3.2g.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-48.253, -45.643, -42.868], [ 41.65 ,  41.087,  45.362]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
    def test_rna_6point4g(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_minmax  = self.o.calc_minmax_all_steps(dcdDataPath+'rna-6.4g.dcd')
	numpy.set_printoptions(precision=3)
	print result_minmax
        expected_minmax = [[-48.253, -45.643, -42.868], [ 41.65 ,  41.332,  45.362]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1],3)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

