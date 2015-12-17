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

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sascalc as sascalc

import numpy

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sascalc_Prop_calcmass(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)
        self.tol = 3

    def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
              self.assertAlmostEqual(a[i],b[i],places)

    def test_null(self):
        with self.assertRaises(Exception):
           self.o.calcmass()

    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [14.00672]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)

    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_mass = [14.00672, 12.01078, 12.01078, 15.99943, 12.01078, 12.01078, 12.01078, \
                          12.01078, 14.00672, 12.01078, 12.01078, 15.99943, 12.01078, 15.99943, 12.01078]
        expected_totalmass = sum(expected_mass)
        self.assertAlmostEqual(expected_totalmass, result_totalmass)
        self.assert_list_almost_equal(expected_mass, result_mass)

    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_totalmass = 106197.087
        self.assertAlmostEqual(expected_totalmass, result_totalmass, self.tol)

    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_totalmass = 4412.904
        self.assertAlmostEqual(expected_totalmass, result_totalmass, self.tol)

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        result_totalmass  = self.o.calcmass()
        result_mass = self.o.mass()
        expected_totalmass = 766109.266
        self.assertAlmostEqual(expected_totalmass, result_totalmass, self.tol)

    def tearDown(self):
        pass 

if __name__ == '__main__': 
   main() 

