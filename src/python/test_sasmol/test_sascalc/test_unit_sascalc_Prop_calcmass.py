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

'''
inf nan and etc need not be tested for this module
'''
from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol

import numpy

class Test_sascalc_Prop_calcmass(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

    def test_null(self):
        self.o.setElement([])
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = []
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)

    def test_one_atoms(self):
        self.o.setElement(['C'])
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [12.01078]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)

    def test_two_atoms(self):
        self.o.setElement(['C','O'])
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [12.01078, 15.99943]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)
       
    def test_six_atoms_duplicate(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [12.01078, 15.99943, 12.01078, 1.00782503214, 83.7982, 107.86822]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)

    def test_wrong_element(self):
        self.o.setElement(['XX','MM'])
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [0.0, 0.0]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)



    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

