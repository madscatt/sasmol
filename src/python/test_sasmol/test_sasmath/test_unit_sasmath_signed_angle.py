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

class Test_sasmath_signed_angle(MockerTestCase): 

    def setUp(self):
        pass 

    def calc_expected(self,a,b,c):
        ada = sum([e*e for e in a])
        bdb = sum([e*e for e in b])

        if numpy.dot(a,a)*numpy.dot(b,b)==0.0:
           return util.num_to_floattype(180.0, floattype)
        else:
           angle = (180.0/util.num_to_floattype(numpy.pi,floattype)) * numpy.arccos(numpy.dot(a,b)/numpy.sqrt(numpy.dot(a,a)*numpy.dot(b,b)))
        sign = cmp(numpy.dot(numpy.cross(a,b),c), 0.0)
        if sign==0: sign=1
        return util.num_to_floattype(sign*angle, floattype)

    def test_all_zero_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 0.0, 0.0],floattype)
        c=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(180.0, floattype)
        self.assertAlmostEqual(result,expected)

    def test_zero_c_array(self):
        a=numpy.array([1.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 1.0, 0.0],floattype)
        c=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(0.0, floattype)
        self.assertAlmostEqual(result,expected)

    def test_unit_arrays(self):
        a=numpy.array([1.0, 0.0, 0.0],floattype)
        b=numpy.array([1.0, 1.0, 0.0],floattype)
        c=numpy.array([0.0, 0.0, 1.0],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(45.0, floattype)
        self.assertAlmostEqual(result,expected)

    def test_same_plane(self):
        a=numpy.array([1.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 1.0, 0.0],floattype)
        c=numpy.array([1.0, 0.0, 0.0],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(0.0, floattype)
        self.assertAlmostEqual(result,expected)

    def test_arb(self):
        a=numpy.array([1.0168308, -1.35572028, -1.35362422],floattype)
        b=numpy.array([-0.69958848, 1.66901076, 0.49978462],floattype)
        #a=numpy.array([1.3, 0.8, -12.8],floattype)
        #b=numpy.array([-9.76, 12.2, 3.0],floattype)
        c=numpy.array([0.3, 0.0, -1.2],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = self.calc_expected(a,b,c)
        self.assertAlmostEqual(result,expected)
        print 'numpy.dot(a,a) ',numpy.dot(a,a)
        print 'numpy.dot(b,b) ',numpy.dot(b,b)
        print 'numpy.dot(a,b) ',numpy.dot(a,b)


    def test_inf_1(self):
        a=numpy.array([util.HUGE, 3.0, 0.0],floattype)
        b=numpy.array([22.3, util.HUGE, 2.3],floattype)
        c=numpy.array([util.HUGE, 1.0, 29.02],floattype)
        result = sasmath.signed_angle(a,b,c)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_inf_2(self):
        a=numpy.array([util.INF, 3.0, 0.0],floattype)
        b=numpy.array([util.INF, 1.2, 2.3],floattype)
        c=numpy.array([0.3, util.INF, 29.02],floattype)
        result = sasmath.signed_angle(a,b,c)
        print result
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_nan(self):
        a=numpy.array([util.NAN, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.NAN, 2.0],floattype)
        c=numpy.array([util.NAN, 1.0, 0.0],floattype)
        result = sasmath.signed_angle(a,b,c)
        self.assertTrue(numpy.isnan(result))

    def test_nan_inf(self):
        a=numpy.array([util.NAN, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.INF, util.NAN],floattype)
        c=numpy.array([util.NAN, 1.0, util.INF],floattype)
        result = sasmath.signed_angle(a,b,c)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_tiny(self):
        a=numpy.array([util.TINY, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.TINY, 2.0],floattype)
        c=numpy.array([2.0, 3.2, util.TINY],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(90.0, floattype)
        self.assertAlmostEqual(result,expected)

    def test_zero(self):
        a=numpy.array([util.ZERO, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.ZERO, 2.0],floattype)
        c=numpy.array([6.7, 1.0, util.ZERO],floattype)
        result = sasmath.signed_angle(a,b,c)
        expected = util.num_to_floattype(90.0, floattype)
        self.assertAlmostEqual(result,expected)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

