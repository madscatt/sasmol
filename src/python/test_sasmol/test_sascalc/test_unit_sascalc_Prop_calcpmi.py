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
import sasmol.sasop as sasop

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_sascalc_Prop_calcpmi(MockerTestCase): 

    def setUp(self):
        self.centertmp = sasop.Move.center

        self.m = Mocker()
        sasop.Move.center = self.m.mock()
        sasop.Move.center(ARGS)
        self.m.result(None)
        self.m.count(0,None)

        self.m.replay()

        self.o=sasmol.SasMol(0)

    def assert_list_almost_equal_flip_sign_allowed(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           sign=1
           for i in range(len(a)):
              if isinstance(a[i],(int,float)):
                 if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
                 if (a[i]*b[i]<0.0): sign = -1
                 self.assertAlmostEqual(a[i],sign*b[i],places)
              else:
                 self.assert_list_almost_equal_flip_sign_allowed(a[i],b[i],places)

    def reorder_eigens(self, result_eigenvalues, result_eigenvectors):
        idx=result_eigenvalues.argsort()
        idx=idx[::-1]
        result_eigenvalues = result_eigenvalues[idx]
        result_eigenvectors = result_eigenvectors[idx]
        result_eigenvectors[2]*=-1
        return result_eigenvalues, result_eigenvectors


    def test_one_atom(self):
        return
        '''
        
        self.o.setCoor(numpy.array([[[-1.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([[156.14, 24.022, 36.032], [24.022, 120.108, -72.065], [36.032, -72.065, 60.054]], floattype)
        expected_eigenvalues = numpy.array([168.151, 168.151, -5.329e-15], floattype)
        expected_eigenvectors = numpy.array([[0.103, -0.812, 0.575], [0.964, 0.148, 0.222], [-0.267, 0.535, 0.802]], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 3)        
        #self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        #self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)
        '''

    def test_two_centered_atoms(self):
        return
        '''
        self.o.setCoor(numpy.array([[[-1.0, -2.0, -3.0],[1.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C','C'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([[26.,  -4.,  -6.], [-4.,  20., -12.], [-6., -12.,  10.]], floattype)     
        expected_eigenvalues = numpy.array([336.302, 336.302, -7.105e-15], floattype)
        expected_eigenvectors = numpy.array([[-0.103, -0.812, 0.575], [0.964, -0.148, -0.222], [0.267, 0.535, 0.802]],floattype)
        #self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        #self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

        '''


    def test_two_uncentered_atoms(self):
        self.o.setCoor(numpy.array([[[-2.0, -2.0, -3.0],[1.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C','N'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print result_I, '\n', result_eigenvalues, '\n',result_eigenvectors
        expected_eigenvalues = numpy.array([400.277, 394.737, 5.54], floattype)
        expected_eigenvectors = numpy.array([[-6.274e-15, -8.321e-01, 5.547e-01], [9.246e-01, -2.114e-01, -3.170e-01], [3.810e-01, 5.129e-01, 7.693e-01]], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

    def test_six_uncentered_atoms(self):
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C','N','O','C','N','O'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print result_I, '\n',result_eigenvalues, '\n',  result_eigenvectors
        expected_eigenvalues = numpy.array([5761.418, 5625.53, 139.66], floattype)
        expected_eigenvectors = numpy.array([[0.351, -0.821, 0.451], [-0.837, -0.059, 0.544],[0.42, 0.568, 0.708]],floattype);
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,2)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)


    def test_six_uncentered_atoms_inf1(self):
        self.o.setCoor(numpy.array([[[util.HUGE, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setMass([1.0, 2.0, 3.2, 3.6, 5.2, 2.8])
        self.o.setNatoms(len(self.o.mass()))
        with self.assertRaises(Exception):
            result = self.o.calcpmi(0)


    def test_six_uncentered_atoms_inf2(self):
        self.o.setCoor(numpy.array([[[util.INF, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setMass([1.0, 2.0, 3.2, 3.6, 5.2, 2.8])
        self.o.setNatoms(len(self.o.mass()))
        with self.assertRaises(Exception):
            result = self.o.calcpmi(0)


    def test_six_uncentered_atoms_nan(self):
        self.o.setCoor(numpy.array([[[util.NAN, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setMass([1.0, 2.0, 3.2, 3.6, 5.2, 2.8])
        self.o.setNatoms(len(self.o.mass()))
        with self.assertRaises(Exception):
            result = self.o.calcpmi(0)


    def test_six_uncentered_atoms_tiny(self):
        self.o.setCoor(numpy.array([[[util.TINY, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C','N','O','C','N','O'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n', list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([[4675.176, -1324.189, -1572.26 ], [-1324.189,  3932.916, -2256.545], [-1572.26 , -2256.545,  2894.494]], floattype)
        expected_eigenvalues = numpy.array([5748.699, 5591.441, 162.447], floattype)
        expected_eigenvectors = numpy.array([[0.321, -0.821, 0.472], [-0.852, -0.032, 0.523], [0.414, 0.57, 0.709]], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,2)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)


    def test_six_uncentered_atoms_ZERO(self):
        self.o.setCoor(numpy.array([[[util.ZERO, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setElement(['C','N','O','C','N','O'])
        self.o.setNatoms(len(self.o.element()))
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([[4675.176, -1324.189, -1572.26 ], [-1324.189,  3932.916, -2256.545], [-1572.26 , -2256.545,  2894.494]], floattype)
        expected_eigenvalues = numpy.array([5748.699, 5591.441, 162.447], floattype)
        expected_eigenvectors = numpy.array([[ 0.321, -0.821,  0.472], [-0.852, -0.032,  0.523], [ 0.414,  0.57 ,  0.709]], floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,2)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

    def tearDown(self):
        self.m.restore()
        self.m.verify()
        
        sasop.Move.center = self.centertmp

if __name__ == '__main__': 
   main() 

