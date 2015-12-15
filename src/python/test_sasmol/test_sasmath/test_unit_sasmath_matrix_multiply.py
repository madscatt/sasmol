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
from mocker import Mocker, MockerTestCase
import sasmol.sasmath as sasmath

import numpy

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_sasmath_matrix_multiply(MockerTestCase): 

    def setUp(self):
        pass 

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

    def test_all_zero_arrays_1(self):
        a=numpy.array([[0.0, 0.0, 0.0]],floattype)
        b=numpy.array([[0.0, 0.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[0.0]],floattype)
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_all_zero_arrays_2(self):
        a=numpy.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]],floattype)
        b=numpy.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[0.0, 0.0],[0.0, 0.0]],floattype)
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_all_zero_arrays_3(self):
        a=numpy.array([[0.0, 0.0, 0.0],[0.0, 0.0, 0.0]],floattype)
        b=numpy.array([[0.0, 0.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[0.0],[0.0]],floattype)
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_unit_arrays_1(self):
        a=numpy.array([[1.0, 0.0, 0.0]],floattype)
        b=numpy.array([[1.0, 1.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[1.0]],floattype)
        self.assert_list_almost_equal(result,expected)
    	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_unit_arrays_2(self):
        a=numpy.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]],floattype)
        b=numpy.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[1.0, 0.0],[0.0, 1.0]],floattype)
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_unit_arrays_3(self):
        a=numpy.array([[1.0, 0.0, 0.0]],floattype)
        b=numpy.array([[1.0, 0.0, 0.0],[0.0, 1.0, 0.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[1.0, 0.0]],floattype)
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_arb_1(self):
        a=numpy.array([[12.0, -20.0, -80.0],[2.02, -901.0, 0.0]],floattype)
        b=numpy.array([[1.23, 0.03, 20.0],[10.0, 1.0, 3.0]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[-1585.84, -140.0],[-24.5454, -880.8]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_arb_2(self):
        a=numpy.array([[1.0168308, -1.35572028, -1.35362422],[-0.69958848, 1.66901076, 0.49978462]],floattype)
        b=numpy.array([-20.69958848, 16.66901076, 20.49978462],floattype)
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[-71.396],[52.547]],floattype)
	print result
        self.assert_list_almost_equal(result,expected,3)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_inf_1(self):
        a=numpy.array([[util.HUGE, -20.0, 80.0],[2.02, util.HUGE, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.HUGE]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[util.INF, util.INF],[util.INF, util.INF]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_inf_2(self):
        a=numpy.array([[util.INF, -20.0, 80.0],[2.02, util.INF, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.INF]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[util.INF, util.INF],[util.INF, util.INF]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_nan(self):
        a=numpy.array([[util.NAN, -20.0, 80.0],[2.02, util.NAN, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.NAN]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[util.NAN, util.NAN],[util.NAN, util.NAN]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_nan_inf(self):
        a=numpy.array([[util.INF, -20.0, 80.0],[2.02, util.NAN, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.NAN]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[util.INF, util.NAN],[util.NAN, util.NAN]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_tiny(self):
        a=numpy.array([[util.TINY, -20.0, 80.0],[2.02, util.TINY, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.TINY]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[1560.0, -420.0],[402.4846, 20.2]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)

    def test_zero(self):
        a=numpy.array([[util.ZERO, -20.0, 80.0],[2.02, util.ZERO, 20.0]],floattype)
        b=numpy.array([[1.23, 2.0, 20.0],[10.0, 21.0, util.ZERO]],floattype).T
        result = sasmath.matrix_multiply(a,b)[1]
        expected = numpy.array([[1560.0, -420.0],[402.4846, 20.2]],floattype)
	print result
        self.assert_list_almost_equal(result,expected)
	result_error = sasmath.matrix_multiply(a,b)[0]
	expected_error = []
	self.assertEqual(result_error, expected_error)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

