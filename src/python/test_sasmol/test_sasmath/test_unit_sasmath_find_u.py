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
import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

PdbPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
modulePdbPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasmath')+os.path.sep

class Test_sasmath(MockerTestCase): 

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

    def assert_list_almost_equal(self,a,b):
        if (len(a)!=len(b)):
           raise "LengthError"
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places=2)

    def test_against_mathematica_two_unit_atoms_one_origin(self):
        x=numpy.array([[0.0, 0.0, 0.0],[1.0, 0.0, 0.0]], floattype)
        y=numpy.array([[0.0, 0.0, 0.0],[0.0, 1.0, 0.0]], floattype)
        result_u = sasmath.find_u(x,y)
        expected_u = [[0.0, 1.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]];
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_against_mathematica_two_unit_atoms(self):
        x=numpy.array([[1.0, 1.0, 1.0], [2.0, 1.0, 1.0]], floattype)
        y=numpy.array([[1.0, 1.0, 1.0], [1.0, 2.0, 1.0]], floattype)
        result_u = sasmath.find_u(x,y); print result_u
        expected_u = [[0.0, 1.0, 0.0], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]];
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_against_mathematica_two_overlap_unit_atoms(self):
        x=numpy.array([[1.0, 1.0, 1.0], [1.0, 2.0, 1.0]], floattype)
        y=numpy.array([[1.0, 1.0, 1.0], [1.0, 2.0, 1.0]], floattype)
        result_u = sasmath.find_u(x,y); print result_u
        print result_u
        expected_u = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_against_mathematica_two_overlap_atoms(self):
        x=numpy.array([[2.920, -2.367, 1.693], [-0.770, -0.827, -0.417]], floattype)
        y=numpy.array([[2.920, -2.367, 1.693], [-0.770, -0.827, -0.417]], floattype)
        result_u = sasmath.find_u(x,y)
        expected_u = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_against_mathematica_three_arbitary_atoms(self):
        x=numpy.array([[2.920, -2.367, 1.693], [-0.770, -0.827, -0.417], [-2.150, 3.193, -1.277]], floattype)
        y=numpy.array([[1.663, -1.170, 3.567], [-1.197, -1.460, -0.523], [-0.467, 2.630, -3.043]], floattype)
        result_u = sasmath.find_u(x,y)
        expected_u = [[0.902737, -0.0539463, 0.426797], [0.23798, 0.8891, -0.390981], [-0.358373, 0.454522, 0.815462]]
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_against_mathematica_six_arbitary_atoms(self):
        x=numpy.array([[0.357, -12.123, 2.098], [1.209, 10.209, -50.082], [-1.098, 3.572, 2.982], \
                       [1.231, -1.230, 0.589], [12.398, -30.289, 19.482], [12.123, 0.980, 19.309]], floattype)
        y=numpy.array([[90.380, 12.987, 0.392], [3.219, 83.390, 0.028], [0.002, 10.298, -18.820], \
                       [12.879, -10.298, 0.987], [0.986, 12.984, 0.367], [12.359, -12.402, 1.298]], floattype)
        result_u = sasmath.find_u(x,y)
        expected_u = [[0.121253, 0.025345, 0.992298], [-0.992602, 0.00937959, 0.12105], [-0.00623933, -0.999635, 0.0262948]]
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_find_u_zero(self):
        x=numpy.array([[0,0,0],[0,0,0]], floattype)
        y=numpy.array([[0,0,0],[0,0,0]], floattype)
        rmx = sasmath.find_u(x,y)
        result = numpy.dot(rmx,y.T)
        result = result.T
        self.assertEqual(len(result),len(x))
        for i in range(len(x)):
          self.assert_list_almost_equal(list(result[i]),x[i])

    def test_find_u_unit(self):
        x=numpy.array([[0,0,0],[1,0,0]], floattype)
        y=numpy.array([[0,0,0],[0,1,0]], floattype)
        rmx = sasmath.find_u(x,y)
        result = numpy.dot(rmx,y.T)
        result = result.T
        self.assertEqual(len(result),len(x))
        for i in range(len(x)):
          self.assert_list_almost_equal(list(result[i]),x[i])

    def test_find_u_arb(self):
        x=numpy.array([[2.920, -2.367, 1.693], [-0.770, -0.827, -0.417], [-2.150, 3.193, -1.277]], floattype)
        y=numpy.array([[1.663, -1.170, 3.567], [-1.197, -1.460, -0.523], [-0.467, 2.630, -3.043]], floattype)
        result_u = sasmath.find_u(x,y)
        expected_u = [[0.902737, -0.0539463, 0.426797], [0.23798, 0.8891, -0.390981], [-0.358373, 0.454522, 0.815462]]
        for i in range(len(result_u)):
          self.assert_list_almost_equal(list(result_u[i]),expected_u[i])

    def test_find_u_rotate_pdb(self):
        m1 = sasmol.SasMol(0)
        m2 = sasmol.SasMol(1)
        m1.read_pdb(PdbPath+"1CRN.pdb")
        m2.read_pdb(modulePdbPath+"1CRN-rot.pdb")

        coor_sub_m1 = m1.coor()[0]
        coor_sub_m2 = m2.coor()[0]

        u = sasmath.find_u(coor_sub_m1, coor_sub_m2)
        result = numpy.array((numpy.matrix(u)*(numpy.matrix(coor_sub_m2).T)).T, numpy.float)
        #print numpy.dot(result.reshape(1,-1)[0],coor_sub_m1.reshape(1,-1)[0])/numpy.sqrt(numpy.dot(coor_sub_m1.reshape(1,-1)[0],coor_sub_m1.reshape(1,-1)[0])*numpy.dot(result.reshape(1,-1)[0],result.reshape(1,-1)[0]))
        #m3 = sasmol.SasMol(2)
        #m3.read_pdb('1CRN.pdb')
        #m3._coor[0,:]=result
        #m3.writepdb('1CRN-result.pdb',0,'w')
        self.assertEqual(len(result),len(coor_sub_m1))
        for i in range(len(coor_sub_m1)):
           self.assert_list_almost_equal(list(result[i]),coor_sub_m1[i])




    def tearDown(self):
        self.m.verify()

if __name__ == '__main__': 
   main() 

