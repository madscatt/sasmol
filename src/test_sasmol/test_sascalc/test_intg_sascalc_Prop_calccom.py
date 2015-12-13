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
from sassie.sasmol import sascalc, sasmol

import numpy

import os

DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'

class Test_sascalc_Prop_calccom(MockerTestCase): 

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
          self.o.calccom(0)


    def test_one_atom_pdb(self):
        self.o.read_pdb(DataPath+'1ATM.pdb')
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [73.944, 41.799, 41.652]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_two_aa_pdb(self):
        self.o.read_pdb(DataPath+'2AAD.pdb')
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [75.68045, 43.70790, 41.27621]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_rna_pdb(self):
        self.o.read_pdb(DataPath+'rna.pdb')
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [-8.033, 4.352, 9.231]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def test_1CRN_pdb(self):
        self.o.read_pdb(DataPath+'1CRN.pdb')
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        expected_com = [9.30026, 9.77488, 6.97776]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_pdb(self):
        self.o.read_pdb(DataPath+'1KP8.pdb')
        self.o.setTotalmass(0.0)
        result_com  = self.o.calccom(0)
        print result_com
        expected_com = [83.286,  0.251, 26.234]
        self.assert_list_almost_equal(expected_com, result_com, self.tol)

    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

