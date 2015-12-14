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
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasproperties')+os.path.sep

class Test_unit_sasio_Files_get_resnames(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      (self.protein, self.dna, self.rna, self.nucleic, self.water)=self.o.get_resnames()

   def unique(self,seq):
      '''
      remove the duplicate elements
      '''
      seen = set()
      seen_add = seen.add
      return [ x for x in seq if x not in seen and not seen_add(x)]

   def compare_namelists(self, l1, l2):
      '''
      check whether two list contains the same atoms
      the order doesnt matter
      '''
      if len(l1)!=len(l2):
         return -1
      else:
         for i1 in l1:
            if i1 not in l2:
               return -1
      return 0


   def test_uniqe(self):
      '''
	   make sure the list is unique
      '''
      #
      unique_protein = self.unique(self.protein)
      self.assertEqual(unique_protein, self.protein)
      unique_dna = self.unique(self.dna)
      self.assertEqual(unique_dna, self.dna)
      unique_rna = self.unique(self.rna)
      self.assertEqual(unique_rna, self.rna)
      unique_nucleic = self.unique(self.nucleic)
      self.assertEqual(unique_nucleic, self.nucleic)
      unique_water = self.unique(self.water)
      self.assertEqual(unique_water, self.water)



   def test_protein(self):
      '''
      make sure the right protein list was generated
      '''
      #
      protein_resnames=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HSD','HSE','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
      self.assertEqual(self.compare_namelists(protein_resnames, self.protein), 0)


   def test_dna(self):
      '''
      make sure the right dna list was generated
      '''
      #
      dna_resnames=['NUSA','NUSG','NUSC','NUSU','DA','DG','DC','DT','ADE','GUA','CYT','THY']
      self.assertEqual(self.compare_namelists(dna_resnames, self.dna), 0)


   def test_rna(self):
      '''
      make sure the right rna list was generated
      '''
      #
      rna_resnames=['RNUS','RNUA','RUUG','RNUC','A', 'C', 'G', 'U','ADE','CYT','GUA','URA']
      self.assertEqual(self.compare_namelists(rna_resnames, self.rna), 0)


   def test_nucleic(self):
      '''
      make sure the right nucleic list was generated
      '''
      #
      nucleic_resnames = ['GUA','ADE','CYT','THY','URA','G', 'A', 'C', 'T', 'U','DA','DG','DC','DT']
      self.assertEqual(self.compare_namelists(nucleic_resnames, self.nucleic), 0)


   def test_water(self):
      '''
      make sure the right water list was generated
      '''
      #
      water_resnames=['TIP3','SPCE','TIP','SPC','TIP4','TP3M'] 
      self.assertEqual(self.compare_namelists(water_resnames, self.water), 0)



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

