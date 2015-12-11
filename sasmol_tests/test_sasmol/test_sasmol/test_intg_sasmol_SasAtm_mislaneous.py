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

'''
sasio.Files.read_pdb seems not getting the moltype right
'''
from unittest import main 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
from sassie.sasmol import sasmol, sasop, sascalc

import numpy, os, copy


import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASSIE_FLOATTYPE']

import os; DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/pdb_common/'


class Test_intg_sasmol_SasAtm_Type(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasAtm(3,'1CRN-3frames.pdb')

   def test_energy(self):
      expected = 10.0
      self.o.setEnergy(expected)
      result = self.o.energy()
      self.assertEqual(expected, result)

   def test_formula(self):
      expected = 'fom'
      self.o.setFormula(expected)
      result = self.o.formula()
      self.assertEqual(expected, result)

   def test_mass(self):
      expected = 100.0
      self.o.setMass(expected)
      result = self.o.mass()
      self.assertEqual(expected, result)

   def test_totalmass(self):
      expected = 100.0
      self.o.read_pdb(DataPath+'1ATM.pdb')
      self.o.setTotalmass(expected)
      result = self.o.totalmass()
      self.assertEqual(expected, result)

   def test_unitcell(self):
      expected = 'P1'
      self.o.setUnitcell(expected)
      result = self.o.unitcell()
      self.assertEqual(expected, result)

   def test_com(self):
      expected = [1.,2.,3.]
      self.o.setCom(expected)
      result = self.o.com()
      self.assertEqual(expected, result)

   def test_natoms(self):
      expected = 100
      self.o.setNatoms(expected)
      result = self.o.natoms()
      self.assertEqual(expected, result)

   def test_rg(self):
      expected = 100.0
      self.o.setRg(expected)
      result = self.o.rg()
      self.assertEqual(expected, result)

   def test_pmi(self):
      expected = 100.0
      self.o.setPmi(expected)
      result = self.o.pmi()
      self.assertEqual(expected, result)

   def test_minimum(self):
      expected = 1.0
      self.o.setMinimum(expected)
      result = self.o.minimum()
      self.assertEqual(expected, result)

   def test_maximum(self):
      expected = 10000.0
      self.o.setMaximum(expected)
      result = self.o.maximum()
      self.assertEqual(expected, result)

   def test_shape(self):
      expected = 'sphere'
      self.o.setShape(expected)
      result = self.o.shape()
      self.assertEqual(expected, result)

   def test_moltype(self):
      expected = 'protein'
      self.o.setMoltype(expected)
      result = self.o.moltype()
      self.assertEqual(expected, result)

   def test_number_of_names(self):
      expected = 100
      self.o.setNumber_of_names(expected)
      result = self.o.number_of_names()
      self.assertEqual(expected, result)

   def test_number_of_resnames(self):
      expected = 100
      self.o.setNumber_of_resnames(expected)
      result = self.o.number_of_resnames()
      self.assertEqual(expected, result)

   def test_number_of_resids(self):
      expected = 100
      self.o.setNumber_of_resids(expected)
      result = self.o.number_of_resids()
      self.assertEqual(expected, result)

   def test_number_of_chains(self):
      expected = 100
      self.o.setNumber_of_chains(expected)
      result = self.o.number_of_chains()
      self.assertEqual(expected, result)

   def test_number_of_segnames(self):
      expected = 100
      self.o.setNumber_of_segnames(expected)
      result = self.o.number_of_segnames()
      self.assertEqual(expected, result)

   def test_number_of_occupancies(self):
      expected = 100
      self.o.setNumber_of_occupancies(expected)
      result = self.o.number_of_occupancies()
      self.assertEqual(expected, result)

   def test_number_of_betas(self):
      expected = 100
      self.o.setNumber_of_betas(expected)
      result = self.o.number_of_betas()
      self.assertEqual(expected, result)

   def test_number_of_elements(self):
      expected = 100
      self.o.setNumber_of_elements(expected)
      result = self.o.number_of_elements()
      self.assertEqual(expected, result)

   def test_names(self):
      expected = ['C','N']
      self.o.setNames(expected)
      result = self.o.names()
      self.assertEqual(expected, result)

   def test_resnames(self):
      expected = ['Ala','Phe']
      self.o.setResnames(expected)
      result = self.o.resnames()
      self.assertEqual(expected, result)

   def test_resids(self):
      expected = [1,2,3]
      self.o.setResids(expected)
      result = self.o.resids()
      self.assertEqual(expected, result)

   def test_chains(self):
      expected = ['A','B']
      self.o.setChains(expected)
      result = self.o.chains()
      self.assertEqual(expected, result)

   def test_segnames(self):
      expected = ['A','A']
      self.o.setSegnames(expected)
      result = self.o.segnames()
      self.assertEqual(expected, result)

   def test_occupancies(self):
      expected = ['1.0','1.0']
      self.o.setOccupancies(expected)
      result = self.o.occupancies()
      self.assertEqual(expected, result)

   def test_betas(self):
      expected = [10.0,10.0]
      self.o.setBetas(expected)
      result = self.o.betas()
      self.assertEqual(expected, result)

   def test_elements(self):
      expected = ['C','N']
      self.o.setElements(expected)
      result = self.o.elements()
      self.assertEqual(expected, result)

   def test_names_mask(self):
      expected = [1,0,1]
      self.o.setNames_mask(expected)
      result = self.o.names_mask()
      self.assertEqual(expected, result)

   def test_resnames_mask(self):
      expected = [1,0]
      self.o.setResnames_mask(expected)
      result = self.o.resnames_mask()
      self.assertEqual(expected, result)

   def test_resids_mask(self):
      expected = [1,0]
      self.o.setResids_mask(expected)
      result = self.o.resids_mask()
      self.assertEqual(expected, result)

   def test_chains_mask(self):
      expected = [1,0]
      self.o.setChains_mask(expected)
      result = self.o.chains_mask()
      self.assertEqual(expected, result)

   def test_occupanies_mask(self):
      expected = [1,0]
      self.o.setOccupancies_mask(expected)
      result = self.o.occupancies_mask()
      self.assertEqual(expected, result)

   def test_betas_mask(self):
      expected = [1,0]
      self.o.setBetas_mask(expected)
      result = self.o.betas_mask()
      self.assertEqual(expected, result)

   def test_elements_mask(self):
      expected = [1,0]
      self.o.setElements_mask(expected)
      result = self.o.elements_mask()
      self.assertEqual(expected, result)

   def test_segnames_mask(self):
      expected = [1,0]
      self.o.setSegnames_mask(expected)
      result = self.o.segnames_mask()
      self.assertEqual(expected, result)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

