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

from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

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



    def test_1_atom(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_2_atoms_duplicate(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[1.0, 2.0, 3.0], [1.0, 2.0, 3.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_2_atoms(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0], [-1.0, 2.0, -6.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[-1.0, 2.0, -6.0], [1.0, 2.0, 3.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[-4.0, 5.0, 6.0],[7.0, 8.0, -9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[-4.0, 2.0, -9.0], [7.0, 8.0, 6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms_inf1(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[-util.HUGE,5.0, 6.0],[7.0, 8.0, -9.0],[1.0, 3.0, 5.0],[2.0, util.HUGE,6.0],[0.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[-util.HUGE,2.0, -9.0], [7.0, util.HUGE,6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms_inf2(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[-util.INF,5.0, 6.0],[7.0, 8.0, -9.0],[1.0, 3.0, 5.0],[2.0, util.INF,6.0],[0.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        expected_minmax = [[-util.INF,2.0, -9.0], [7.0, util.INF,6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms_nan(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.NAN,5.0, 6.0],[7.0, 8.0, -9.0],[1.0, 3.0, 5.0],[2.0, util.NAN,6.0],[0.0, 2.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        print result_minmax
        expected_minmax = [[util.NAN,util.NAN,-9.0], [util.NAN,util.NAN,6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms_tiny(self):
        self.o.setCoor(numpy.array([[[1.0, -2.0, 3.0],[-util.TINY,-5.0, 6.0],[7.0, -8.0, -9.0],[1.0, -3.0, 5.0],[2.0, util.TINY,6.0],[0.0, 0.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        print result_minmax
        expected_minmax = [[-util.TINY,-8.0, -9.0], [7.0, util.TINY,6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def test_6_atoms_zero(self):
        self.o.setCoor(numpy.array([[[1.0, -2.0, 3.0],[-util.ZERO,-5.0, 6.0],[7.0, -8.0, -9.0],[1.0, -3.0, 5.0],[2.0, util.ZERO,6.0],[0.0, 0.0, 3.0]]],floattype))
        result_minmax  = self.o.calcminmax()
        print result_minmax
        expected_minmax = [[-util.ZERO,-8.0, -9.0], [7.0, util.ZERO,6.0]]
        self.assert_list_almost_equal(expected_minmax[0], result_minmax[0])
        self.assert_list_almost_equal(expected_minmax[1], result_minmax[1])

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

