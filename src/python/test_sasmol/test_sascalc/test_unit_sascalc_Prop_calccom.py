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

class Test_sascalc_Prop_calccom(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)
        self.tol = 3

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

    def test_null(self):
        self.o.setElement([])
        self.o.setCoor(numpy.zeros((1.0, 0.0, 3.0),floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [util.NAN, util.NAN, util.NAN]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_one_atom(self):
        self.o.setElement(['C'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [1.0, 2.0, 3.0]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_two_atoms(self):
        self.o.setElement(['C', 'AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [3.69943, 4.69942, 5.69943]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [1.41253, 3.24053, 4.60498]
        print '\nresult_com \n',result_com
        print '\nexpected_com \n',expected_com
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate_inf_1(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[util.HUGE, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [util.INF, 3.24054, 4.60499]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate_inf_2(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.INF, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [util.INF, 3.24054, 4.60499]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate_nan(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, util.NAN, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [1.41253, util.NAN, 4.60499]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate_tiny(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.TINY, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [1.13750, 3.24053, 4.60499]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_six_atoms_duplicate_zero(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.ZERO, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [1.13750, 3.24053, 4.60499]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_wrong_element(self):
        self.o._element = ['X','M']
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0]]],floattype))
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [util.NAN, util.NAN, util.NAN]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

