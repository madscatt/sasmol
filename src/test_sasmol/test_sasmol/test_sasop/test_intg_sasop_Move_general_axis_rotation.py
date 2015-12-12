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



from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS
from sassie.sasmol import sasmol, sasop, sascalc

import numpy


import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'

class Test_intg_sasop_general_rotate(MockerTestCase): 

    def setUp(self):
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


    def test_one_atom_pdb_norotate(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        axis = 'x'
        frame = 0
        theta=0.0
        #
        self.o.general_axis_rotate(frame,theta,0.2,1.3,-3.5)
        result_coor = self.o.coor()
        result_com  = self.o.calccom(0)
        print '\nresult_coor:\n'; util.printfl([result_coor]); print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_coor = numpy.array([[[73.944, 41.799, 41.652]]], floattype)
        expected_com = numpy.array([73.944, 41.799, 41.652], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)


    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        axis = 'x'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.general_axis_rotate(frame,theta,0.2,1.3,-3.5)
        result_coor = self.o.coor()
        result_com  = self.o.calccom(0)
        print '\nresult_coor:\n'; util.printfl([result_coor]); print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_coor = numpy.array([[[-215.775, 167.484, 356.058]]], floattype)
        expected_com = numpy.array([-215.775, 167.484, 356.058], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)


    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.general_axis_rotate(frame,theta,0.2,1.3,-3.5)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([-221.139, 178.873, 343.429], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        frame = 0
        theta=numpy.pi/3.0
        #
        self.o.general_axis_rotate(frame,theta,0.2,1.3,-3.5)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([-30.425, -38.942, 44.267], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.general_axis_rotate(frame,theta,0,1,0)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([-6.978, 9.775, 9.300], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        frame = 0
        theta=12.0
        #
        self.o.general_axis_rotate(frame,theta,0,1,0)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([84.358, 0.251, -22.552], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

