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

class Test_sasmath_dihedral_angle(MockerTestCase): 

    def setUp(self):
        pass 



    def test_all_zero_arrays(self):
        a=numpy.array([0.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 0.0, 0.0],floattype)
        c=numpy.array([0.0, 0.0, 0.0],floattype)
        d=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 180.0
        self.assertAlmostEqual(result,expected)

    def test_unit_arrays(self):
        a=numpy.array([1.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 0.0, 0.0],floattype)
        c=numpy.array([0.0, 1.0, 0.0],floattype)
        d=numpy.array([0.0, 0.0, 1.0],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = -90.0
        self.assertAlmostEqual(result,expected)

    def test_same_plane_1(self):
        a=numpy.array([1.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 0.0, 0.0],floattype)
        c=numpy.array([0.0, 1.0, 0.0],floattype)
        d=numpy.array([1.0, 1.0, 0.0],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 0.0
        self.assertAlmostEqual(result,expected)

    def test_arb_1(self):
        a=numpy.array([23.0, 0.0, 0.0],floattype)
        b=numpy.array([0.0, 0.0, 0.0],floattype)
        c=numpy.array([0.0, 12.0, 0.0],floattype)
        d=numpy.array([0.0, 212.0, -20.0],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 90.0
        self.assertAlmostEqual(result,expected)

    def test_arb_2(self):
        a=numpy.array([1.0168308, -1.35572028, -1.35362422],floattype)
        b=numpy.array([-0.69958848, 1.66901076, 0.49978462],floattype)
        c=numpy.array([12.0168308, 12.35572028, 1.35362422],floattype)
        d=numpy.array([-20.69958848, 16.66901076, 20.49978462],floattype)
        c=numpy.array([0.3, 0.0, -1.2],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 85.786
        self.assertAlmostEqual(result,expected,3)

    def test_inf_1(self):
        a=numpy.array([util.HUGE, 3.0, 0.0],floattype)
        b=numpy.array([22.3, util.HUGE, 2.3],floattype)
        c=numpy.array([util.HUGE, 1.0, 29.02],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_inf_2(self):
        a=numpy.array([util.INF, 3.0, 0.0],floattype)
        b=numpy.array([util.INF, 1.2, 2.3],floattype)
        c=numpy.array([0.3, util.INF, 29.02],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_nan(self):
        a=numpy.array([util.NAN, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.NAN, 2.0],floattype)
        c=numpy.array([util.NAN, 1.0, 0.0],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        self.assertTrue(numpy.isnan(result))

    def test_nan_inf(self):
        a=numpy.array([util.NAN, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.INF, util.NAN],floattype)
        c=numpy.array([util.NAN, 1.0, util.INF],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        self.assertTrue(numpy.isnan(result) or numpy.isinf(result))

    def test_tiny(self):
        a=numpy.array([util.TINY, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.TINY, 2.0],floattype)
        c=numpy.array([2.0, 3.2, util.TINY],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 66.032
        self.assertAlmostEqual(result,expected,3)

    def test_zero(self):
        a=numpy.array([util.ZERO, 3.0, 0.0],floattype)
        b=numpy.array([6.0, util.ZERO, 2.0],floattype)
        c=numpy.array([6.7, 1.0, util.ZERO],floattype)
        d=numpy.array([2.1, 1.0, 29.02],floattype)
        result = sasmath.dihedral_angle(a,b,c,d)
        expected = 89.502
        self.assertAlmostEqual(result,expected,3)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

