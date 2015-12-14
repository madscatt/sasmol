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

from sassie.core_testing.util import env, util

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase
from sassie.sasmol import sasmol

import numpy

import os
floattype=os.environ['SASSIE_FLOATTYPE']


PdbPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'

class Test_sascalc_Prop_calcpmi(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

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


    def test_null(self):
        with self.assertRaises(Exception):
            self.o.read_pdb(PdbPath+'NULL.pdb')
        with self.assertRaises(Exception):
           result_pmi  = self.o.calcpmi(0)

    def test_one_atom_pdb(self):
        self.o.read_pdb(PdbPath+'1ATM.pdb')
        self.o.calcmass()
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.zeros((3,3),floattype)
        expected_eigenvalues = numpy.zeros(3,floattype)
        expected_eigenvectors = numpy.array([[0.0, 0.0, 1.0], [0.0, 1.0, 0.0], [-1.0, 0.0, 0.0]],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 5)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)


    def test_two_aa_pdb(self):
        self.o.read_pdb(PdbPath+'2AAD.pdb')
        self.o.calcmass()
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([numpy.array([ 589.53374631,  -64.32846157, -439.3753857 ]), numpy.array([  -64.32846157,  1532.13560848,   -65.3989943 ]), numpy.array([ -439.3753857 ,   -65.3989943 ,  1407.58300946])], floattype)
        expected_eigenvalues = numpy.array([1614.281458830048, 1523.0603348786992, 391.91057054356725], floattype)
        expected_eigenvectors = numpy.array([numpy.array([ 0.33751536,  0.41020629, -0.84723915]), numpy.array([ 0.22717091, -0.90894656, -0.3495848 ]), numpy.array([-0.913497  , -0.07447786, -0.39997036])],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 5)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)


    def test_rna_pdb(self):
        self.o.read_pdb(PdbPath+'rna.pdb')
        self.o.calcmass()
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([numpy.array([  3.04411898e+08,   4.04333713e+07,   4.06520707e+07]), numpy.array([  4.04333713e+07,   3.08529104e+08,  -4.59336765e+07]), numpy.array([  4.06520707e+07,  -4.59336765e+07,   3.02196582e+08])], floattype)
        expected_eigenvalues = numpy.array([351687532.76625204, 343174952.58514869, 220275098.79483908], floattype)
        expected_eigenvectors = numpy.array([numpy.array([-0.1525973 , -0.78373478,  0.60205802]), numpy.array([ 0.8122177 ,  0.24761209,  0.52819567]), numpy.array([ 0.56304216, -0.56960341, -0.59877832])],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, -1)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

    def test_1CRN_pdb(self):
        self.o.read_pdb(PdbPath+'1CRN.pdb')
        self.o.calcmass()
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([numpy.array([ 258450.06343405,  -45258.19954399,  -67627.06498969]), numpy.array([ -45258.19954399,  311061.10234081,   -3089.06334182]), numpy.array([ -67627.06498969,   -3089.06334182,  244380.69561839])], floattype)
        expected_eigenvalues = numpy.array([349987.99722910166, 288718.59037127224, 175185.27379287069], floattype)
        expected_eigenvectors = numpy.array([numpy.array([ 0.61882991, -0.68963421, -0.37610398]), numpy.array([-0.37918551, -0.68156978,  0.62584422]), numpy.array([-0.68794469, -0.24467794, -0.68327506])],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, 3)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(PdbPath+'1KP8.pdb')
        self.o.calcmass()
        result = self.o.calcpmi(0)
        result_eigenvalues = result[0]
        result_eigenvectors = result[1].T
        result_I = result[2]
        result_eigenvalues, result_eigenvectors = self.reorder_eigens(result_eigenvalues, result_eigenvectors)
        print list(result_I), '\n',list(result_eigenvalues), '\n', list(result_eigenvectors)
        expected_I = numpy.array([numpy.array([  2.11885718e+09,  -5.07311719e+06,  -6.58159781e+06]), numpy.array([ -5.07311719e+06,   2.11848735e+09,   7.27900160e+06]), numpy.array([ -6.58159781e+06,   7.27900160e+06,   1.90342505e+09])], floattype)
        expected_eigenvalues = numpy.array([2124183018.8505797, 2113597829.3673213, 1902988729.7621682], floattype)
        expected_eigenvectors = numpy.array([numpy.array([ 0.71720897, -0.69544778, -0.04431345]), numpy.array([-0.69622572, -0.71781641, -0.003058  ]), numpy.array([-0.02968224,  0.03304539, -0.999013  ])],floattype)
        self.assert_list_almost_equal_flip_sign_allowed(expected_I, result_I, -2)        
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvalues, result_eigenvalues,3)
        self.assert_list_almost_equal_flip_sign_allowed(expected_eigenvectors, result_eigenvectors,3)


    def test_problem_pdb(self):
        with self.assertRaises(Exception):
           self.o.read_pdb(PdbPath+'1PSI.pdb')
        with self.assertRaises(Exception):
           result_pmi  = self.o.calcpmi(0)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

