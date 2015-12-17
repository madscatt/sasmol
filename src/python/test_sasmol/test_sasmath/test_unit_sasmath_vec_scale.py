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
import sasmol.sasmath as sasmath
import sasmol.sasio as sasio

import numpy

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_sasmath_vec_scale(MockerTestCase): 

    def setUp(self):
        pass


    def assert_list_almost_equal(self,a,b,places=7):
        if (len(a)!=len(b)):
           raise "LengthError"
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

    def test_zero_array(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.num_to_floattype(0,floattype)
        y=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_unit_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.num_to_floattype(1.0,floattype)
        y=numpy.array([1.0, 1.0, 1.0],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_arb_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.num_to_floattype(-35.01,floattype)
        y=numpy.array([-23.1, 0.98, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_inf_1(self):
        a=numpy.array([0.02, 5.2, -10.0],floattype)
        x=util.HUGE
        y=numpy.array([-23.1, 19.8, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_inf_2(self):
        a=numpy.array([util.HUGE, 103.0, -30.0],floattype)
        x=util.INF
        y=numpy.array([-23.1, 0.98, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_nan(self):
        a=numpy.array([20.0, util.NAN, -30.0],floattype)
        x=util.NAN
        y=numpy.array([-23.1, util.NAN, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_nan_inf(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.INF
        y=numpy.array([-23.1, util.NAN, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_tiny(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.TINY
        y=numpy.array([-23.1, util.TINY, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)

    def test_zero(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=util.ZERO
        y=numpy.array([-23.1, util.ZERO, 18.208],floattype)
        result = sasmath.vec_scale(a,x,y)
        expected = util.list_to_floattype(x*y,floattype)
        self.assert_list_almost_equal(result,expected)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

