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
from mocker import Mocker, MockerTestCase, ANY, ARGS
from sassie.sasmol import sasmol

import numpy


import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

PdbPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'
modulePdbPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/sasmol/sascalc/'

class test_sascalc_Prop_calcrmsd(MockerTestCase): 

    def setUp(self):
        self.o1=sasmol.SasMol(0)
        self.o2=sasmol.SasMol(0)
        self.tol = 3

    def test_null(self):
        with self.assertRaises(Exception):
            self.o1.calcrmsd(self.o2)

    def test_one_atom_pdb(self):
        self.o1.read_pdb(PdbPath+'1ATM.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.read_pdb(modulePdbPath+'1ATN.pdb')
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = 189.20635
        self.assertAlmostEqual(expected_rmsd, result_rmsd,self.tol)

    def test_one_overlap_atom_pdb(self):
        self.o1.read_pdb(PdbPath+'1ATM.pdb')                                                            
        self.o2.read_pdb(PdbPath+'1ATM.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = 0.0
        self.assertAlmostEqual(expected_rmsd, result_rmsd, self.tol)

    def test_1CRN_1CRNrot_pdb(self):
        self.o1.read_pdb(PdbPath+'1CRN.pdb')                                                            
        self.o2.read_pdb(modulePdbPath+'1CRN-rot.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)                                              
        expected_rmsd = 29.008
        self.assertAlmostEqual(expected_rmsd, result_rmsd, self.tol) 

    def test_1CRN_1CRNrotshift_pdb(self):
        self.o1.read_pdb(PdbPath+'1CRN.pdb')                                                            
        self.o2.read_pdb(modulePdbPath+'1CRN-rot-shift.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)                                              
        expected_rmsd = 19.831
        self.assertAlmostEqual(expected_rmsd, result_rmsd, self.tol) 

    def test_rna_rna_pdb(self):
        self.o1.read_pdb(PdbPath+'rna.pdb')                                                            
        self.o2.read_pdb(PdbPath+'rna.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)                                              
        expected_rmsd = 0.0
        self.assertAlmostEqual(expected_rmsd, result_rmsd, self.tol) 

    @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
    def test_1KP8_1KP8_pdb(self):
        self.o1.read_pdb(PdbPath+'1KP8.pdb')                                                            
        self.o2.read_pdb(PdbPath+'1KP8.pdb')
        self.o1.setNatoms(len((self.o1.coor())[0]))
        self.o2.setNatoms(len((self.o2.coor())[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)                                              
        expected_rmsd = 0.0
        self.assertAlmostEqual(expected_rmsd, result_rmsd, self.tol) 


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

