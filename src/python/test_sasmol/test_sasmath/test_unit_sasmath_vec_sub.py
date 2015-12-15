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

class Test_sasmath_vec_sub(MockerTestCase): 

    def setUp(self):
        pass 


    def assert_list_almost_equal(self,a,b):
        if (len(a)!=len(b)):
           raise "LengthError"
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i])

    def test_all_zero_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=numpy.array([0.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_unit_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=numpy.array([1.0, 1.0, 1.0],floattype)
        y=numpy.array([1.0, 1.0, 1.0],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_arb_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, 5.6, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_inf_1(self):
        a=numpy.array([util.HUGE, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, util.HUGE, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_inf_2(self):
        a=numpy.array([util.INF, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, util.INF, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_nan(self):
        a=numpy.array([util.NAN, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, util.NAN, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_inf_nan(self):
        a=numpy.array([util.INF, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, util.NAN, 0.2],floattype)
        y=numpy.array([12.01, -109.3, util.HUGE],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_tiny(self):
        a=numpy.array([util.TINY, 0.0, 0.0],floattype)
        x=numpy.array([-2.3, util.TINY, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)

    def test_zero(self):
        a=numpy.array([util.ZERO, 2.0, 3.0],floattype)
        x=numpy.array([-2.3, util.ZERO, 0.2],floattype)
        y=numpy.array([12.01, -109.3, 3.87],floattype)
        result = sasmath.vec_sub(a,x,y)
        expected = numpy.subtract(x,y)
        self.assert_list_almost_equal(result,expected)


    def tearDown(self):
        pass


if __name__ == '__main__': 
   main() 

