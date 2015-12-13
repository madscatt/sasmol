'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
	 Core-Testing: Copyright (C) 2011 Hailiang Zhang, Ph.D.

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

from sassie.core_testing.util import env,util



from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
from sassie.sasmol import sasmol, sasop, sascalc

import numpy


import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_unit_sasop_Move_moveto(MockerTestCase): 

    def setUp(self):
        self.back_masscheck = sasop.Move.masscheck
        self.back_calccom = sascalc.Prop.calccom 

        self.m = Mocker()

        sasop.Move.masscheck = self.m.mock()
        sasop.Move.masscheck(ARGS)
        self.m.result(None)
        self.m.count(0,None)

        sascalc.Prop.calccom = self.m.mock()
        sascalc.Prop.calccom(ARGS)
        self.m.result(None)
        self.m.count(0,None)

        self.m.replay()

        self.o=sasmol.SasMol(0)

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

    def test_one_atom(self):
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.0]]],floattype))
        self.o.setCom(numpy.array([-1.0, 2.0, 3.0],floattype))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()
        print result_coor
        expected_coor = numpy.array([[[1.0, 3.0, 6.0]]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_two_atoms(self):
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.87],[-5.0, 3.2, 6.0]]],floattype))
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor
        expected_coor = numpy.array([[3.0, 2.4, 4.935],[-1.0, 3.6, 7.065]], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_six_atoms(self):
        coorback = numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_inf1(self):
        coorback = numpy.array([[[1.2, 2.0, util.HUGE],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, util.HUGE],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_inf2(self):
        coorback = numpy.array([[[1.2, 2.0, util.INF],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, util.INF],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_nan(self):
        coorback = numpy.array([[[1.2, 2.0, util.NAN],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, util.NAN],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_tiny(self):
        coorback = numpy.array([[[1.2, 2.0, util.TINY],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, util.TINY],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)

    def test_six_atoms_zero(self):
        coorback = numpy.array([[[1.2, 2.0, util.ZERO],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        coor = numpy.array([[[1.2, 2.0, util.ZERO],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype)
        self.o.setCoor(coor)
        self.o.setCom(numpy.average(self.o.coor()[0],axis=0))
        value = numpy.array([1.0, 3.0, 6.0],floattype)
        self.o.moveto(0,value)
        result_coor = self.o.coor()[0]
        print result_coor,'\n\n'
        expected_coor = coorback[0]-self.o.com()+value
        print expected_coor
        self.assert_list_almost_equal(expected_coor, result_coor,3)



    def tearDown(self):
        self.m.verify()
        sascalc.Prop.calccom  =self.back_calccom
        sasop.Move.masscheck  =self.back_masscheck


if __name__ == '__main__': 
   main() 

