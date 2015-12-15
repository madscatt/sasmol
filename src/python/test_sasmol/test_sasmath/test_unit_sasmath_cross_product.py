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

class Test_sasmath_cross_product(MockerTestCase): 

    def setUp(self):
        self.m = Mocker()

        """
        sasmath.Math.__init__ = self.m.mock()
        sasmath.Math.__init__(ARGS)
        self.m.result(None)
        self.m.count(0,None)
        """

        self.m.replay()

        #self.o=sasmath.Math()

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise "LengthError"
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

    def test_two_zero_arrays(self):

        x=numpy.array([0.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 0.0, 0.0],floattype)
        result = sasmath.cross_product(x,y)
        self.assert_list_almost_equal(a=result,b=[0.0, 0.0, 0.0])

    def test_1_zero_array(self):

        x=numpy.array([0.0, 0.0, 0.0],floattype)
        y=numpy.array([3.2, 2.6, -0.3],floattype)
        result=sasmath.cross_product(x,y)
        self.assert_list_almost_equal(result,[0.0, 0.0, 0.0])

    def test_x_axis_cross_y_axis(self):
        x=numpy.array([1.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 1.0, 0.0],floattype)
        result=sasmath.cross_product(x,y)
        self.assert_list_almost_equal(result,[0.0, 0.0, 1.0])

    def test_x_axis_cross_z_axis(self):
        x=numpy.array([1.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 0.0, 1.0],floattype)
        result=sasmath.cross_product(x,y)
        self.assert_list_almost_equal(result,[0.0, -1.0, 0.0])

    def test_2x_axis_cross_3z_axis(self):
        x=numpy.array([2.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 0.0, 3.0],floattype)
        result=sasmath.cross_product(x,y)
        self.assert_list_almost_equal(result,[0.0, -6.0, 0.0])

    def test_arbitary(self):
        x=numpy.array([3.0, 0.0, 0.0],floattype)
        y=numpy.array([0.0, 0.0, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        self.assert_list_almost_equal(result,[0.0, 6.6, 0.0])

    def test_against_numpy(self):
        x=numpy.array([3.0, 3.0, -100.2],floattype)
        y=numpy.array([17.68, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)

    def test_inf_1(self):
        x=numpy.array([util.HUGE, 3.0, -100.2],floattype)
        y=numpy.array([2.0, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        print result
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)

    def test_inf_2(self):
        x=numpy.array([util.INF, 3.0, -100.2],floattype)
        y=numpy.array([1.0, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)

    def test_nan(self):
        x=numpy.array([util.NAN, 3.0, -100.2],floattype)
        y=numpy.array([1.0, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)

    def test_tiny(self):
        x=numpy.array([util.TINY, 3.0, -100.2],floattype)
        y=numpy.array([1.0, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)

    def test_zero(self):
        x=numpy.array([util.ZERO, 3.0, -100.2],floattype)
        y=numpy.array([1.0, 0.9, -2.2],floattype)
        result=sasmath.cross_product(x,y)
        resultnp = numpy.cross(x,y)
        self.assert_list_almost_equal(result,resultnp)


    def tearDown(self):
        self.m.verify()

if __name__ == '__main__': 
   main() 

