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

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop
import sasmol.sascalc as sascalc

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_intg_sasop_Move_rotate(MockerTestCase): 

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


    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        axis = 'x'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()
        result_com  = self.o.calccom(0)
        print '\nresult_coor:\n'; util.printfl([result_coor]); print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_coor = numpy.array([[[73.944, -41.652, 41.799]]], floattype)
        expected_com = numpy.array([73.944, -41.652, 41.799], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,3)
        self.assert_list_almost_equal(expected_com, result_com,3)


    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        axis = 'y'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.rotate(frame,axis,theta)
        result_coor = self.o.coor()
        result_com  = self.o.calccom(0)
        print '\nresult_coor:\n'; util.printfl([result_coor]); print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_coor = numpy.array([[[41.652, 41.799, -73.944], [40.456, 42.563, -74.229], [40.463, 43.093, -75.667], [39.401, 43.279, -76.264], [40.336, 43.734, -73.210], [39.926, 43.168, -71.856], [39.354, 44.782, -73.67], [39.946, 44.177, -70.721], [41.647, 43.330, -76.231], [41.730, 43.852, -77.592], [42.184, 42.820, -78.617], [42.656, 43.169, -79.712], [42.648, 45.097, -77.671], [43.910, 44.816, -77.054], [42.000, 46.273, -76.970]]], floattype)
        expected_com = numpy.array([41.276, 43.708, -75.680], floattype)
        self.assert_list_almost_equal(expected_coor, result_coor,1)
        self.assert_list_almost_equal(expected_com, result_com,2)


    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        axis = 'z'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.rotate(frame,axis,theta)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([-4.352, -8.033, 9.231], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        axis = 'z'
        frame = 0
        theta=numpy.pi/2.0
        #
        self.o.rotate(frame,axis,theta)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([-9.775, 9.300, 6.978], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)


    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        axis = 'x'
        frame = 0
        theta=12.0
        #
        self.o.rotate(frame,axis,theta)
        result_com  = self.o.calccom(0)
        print '\nresult_com:\n',util.printfl([result_com])
        #
        expected_com = numpy.array([83.286, 14.288, 22.003], floattype)
        self.assert_list_almost_equal(expected_com, result_com,2)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

