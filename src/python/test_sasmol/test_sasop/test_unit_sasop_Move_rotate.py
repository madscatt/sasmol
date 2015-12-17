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

from sasmol.test_sasmol.util import env,util

from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop
import sasmol.sascalc as sascalc

import numpy, os, copy

import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASSIE_FLOATTYPE']

class Test_unit_sasop_Move_translate(MockerTestCase): 

    def setUp(self):
        self.m = Mocker()

        self.back_masscheck = sasop.Move.masscheck
        sasop.Move.masscheck = self.m.mock()
        sasop.Move.masscheck(ARGS)
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


    def test_one_atom_about_X_axis(self):
        '''
	test one atom rotating about the z axis
	'''
        self.o.setCoor(numpy.array([[[1.0, 1.0, 1.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(90.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        expected_coor = numpy.array([[[1.0, -1.0, 1.0]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_one_atom_about_Y_axis(self):
        '''
	test one atom rotating about the z axis
	'''
        self.o.setCoor(numpy.array([[[1.0, 1.0, 1.0]]],floattype))
        axis = 'y'
        theta = util.num_to_floattype(90.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        expected_coor = numpy.array([[[1.0, 1.0, -1.0]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_one_atom_about_Z_axis(self):
        '''
	test one atom rotating about the z axis
	'''
        self.o.setCoor(numpy.array([[[1.0, 1.0, 1.0]]],floattype))
        axis = 'z'
        theta = util.num_to_floattype(90.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        expected_coor = numpy.array([[[-1.0, 1.0, 1.0]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_two_atoms_about_X_axis(self):
	'''
	test two atoms rotating about the x axis
	'''
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.87],[-5.0, 3.2, 6.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(90.0*numpy.pi/180.0, floattype)
        frame = 0
        #
	expected_coor = numpy.array([[[-1.0, -3.87, 2.0],[-5.0, -6.0, 3.2]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_two_atoms_about_Y_axis(self):
	#test two atoms rotating about the x axis
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.87],[-5.0, 3.2, 6.0]]],floattype))
        axis = 'y'
        theta = util.num_to_floattype(90.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        expected_coor = numpy.array([[[3.87, 2.0, 1.0],[6.0, 3.2, 5.0]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_two_atoms_about_Z_axis(self):
        '''
	test two atoms rotating about the x axis
        '''
	self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.87],[-5.0, 3.2, 6.0]]],floattype))
        axis = 'z'
        theta = util.num_to_floattype(280.0*numpy.pi/180.0, floattype)
        frame = 0
        #
	#expected_coor = numpy.array([[[-2.0, -1.0, 3.87],[-3.2, -5.0, 6.0]]], floattype)[frame]
	cs=numpy.cos(theta)
	si=numpy.sin(theta)
	expected_coor = numpy.array([[[-1.0*cs-2.0*si, -1.0*si+2.0*cs, 3.87],[-5.0*cs-3.2*si, -5.0*si+3.2*cs, 6.0]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_six_atoms_about_X_axis(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(128.0*numpy.pi/180.0, floattype)
        frame = 0
        #
	cs=numpy.cos(theta)
	si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, 2.0*cs-3.0*si, 2.0*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_six_atoms_about_X_axis_inf1(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(util.HUGE*numpy.pi/180.0, floattype)
        frame = 0
        #
	cs=numpy.cos(theta)
	si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, 2.0*cs-3.0*si, 2.0*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    #the first element becomes nan after 1.2+inf*0(nan)
    def test_six_atoms_about_X_axis_inf2(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(util.INF*numpy.pi/180.0, floattype)
        frame = 0
        #
        cs=numpy.cos(theta)
        si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, 2.0*cs-3.0*si, 2.0*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    #the first row of elements becomes all nan due to operation between the regular number with nan*0
    def test_six_atoms_about_X_axis_nan(self):
        self.o.setCoor(numpy.array([[[1.2, 2.0, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(util.NAN*numpy.pi/180.0, floattype)
        frame = 0
        #
        cs=numpy.cos(theta)
        si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, util.NAN*cs-3.0*si, util.NAN*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_six_atoms_about_X_axis_tiny(self):
        self.o.setCoor(numpy.array([[[1.2, util.TINY, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(128.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        cs=numpy.cos(theta)
        si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, util.TINY*cs-3.0*si, util.TINY*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)


    def test_six_atoms_about_X_axis_zero(self):
        self.o.setCoor(numpy.array([[[1.2, util.ZERO, 3.0],[-2.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        axis = 'x'
        theta = util.num_to_floattype(128.0*numpy.pi/180.0, floattype)
        frame = 0
        #
        cs=numpy.cos(theta)
        si=numpy.sin(theta)
        expected_coor = numpy.array([[[1.2, util.ZERO*cs-3.0*si, util.ZERO*si+3.0*cs],[-2.0, 5.0*cs-6.0*si, 5.0*si+6.0*cs],[7.0, 8.0*cs-9.0*si, 8.0*si+9.0*cs],[1.0, 3.0*cs-5.0*si, 3.0*si+5.0*cs],[2.0, 4.0*cs-6.0*si, 4.0*si+6.0*cs],[0.0, 2.0*cs-3.0*si, 2.0*si+3.0*cs]]],floattype)[frame]
	print 'expected_coor:\n', expected_coor
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()[frame]
        print 'result_coor:\n',result_coor
        #
        self.assert_list_almost_equal(expected_coor, result_coor,3)



    def tearDown(self):
        self.m.verify()
        sasop.Move.masscheck  =self.back_masscheck



if __name__ == '__main__': 
   main() 

