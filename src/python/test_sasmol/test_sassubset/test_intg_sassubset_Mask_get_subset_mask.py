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

from sasmol.test_sasmol.util import env

"""
Integration test for sasio.sassubset.Mask.get_dihedral_subset_mask

contract:

null test by reading a nonexisting pdb, mask none

negative test by providing the wrong basis_filter

1atom pdb, mask none due to empty basis_filter ('')
1atom pdb, mask none due to illegal basis_filter ('abc')
1atom pdb, mask none due to wrong basis_filter ('name[i]=="B"')
1atom pdb, mask all by basis_filter('name[i]=="N"')
1atom pdb, mask all by basis_filter('resid[i]==515')
1atom pdb, mask all by basis_filter('name[i]=="N" and resid[i]==515')

2 residues and 15 atoms, mask none due to empty basis_filter ('')
2 residues and 15 atoms, mask none due to illegal basis_filter ('abc')
2 residues and 15 atoms, mask none due to wrong basis_filter ('name[i]=="B" or resid[i]==12')
2 residues and 15 atoms, mask 8 by atom name selection
2 residues and 15 atoms, mask 8 by resid selection
2 residues and 15 atoms, mask 8 by elaborate selection

rna molecule, mask none due to empty basis_filter
rna molecule, mask none due to wrong basis_filter('moltype[i]=="protein"')
rna molecule, mask all

small protein (crambin with 46 residues), mask none due to wrong basis_filter
small protein (crambin with 46 residues), mask 4
small protein (crambin with 46 residues), mask all

large protein (groel with 526*14 residues), mask none (Skipped as SASSIE_LARGETEST)
large protein (groel with 526*14 residues), mask partial (Skipped as SASSIE_LARGETEST)
large protein (groel with 526*14 residues), mask all (Skipped as SASSIE_LARGETEST)

problemetic pdb (1PSI wih unpaird MODEL/ENDMDL)
"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset

import numpy

import os

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sassubset_Mask_get_subset_mask(MockerTestCase): 
 

   def setUp(self):
      self.o=sasmol.SasMol(0)


   def test_null(self):
      '''
      test for a null atom list
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      basis_filter = 'name[i]=="N"'
      #
      with self.assertRaises(Exception):
         self.o.get_subset_mask(basis_filter)


   def test_negative(self):
      '''
      negative test by providing the wrong basis_filter
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      basis_filter = 'abc'
      #
      with self.assertRaises(Exception):
         self.o.get_subset_mask(basis_filter)


   def test_1ATM_0outof1_1(self):
      '''
	   test a pdb file with 1 atom
      0/1 will be selected due to empty basis_filter
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = ''
      #
      expecting_error = True
      expectd_mask = []
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1ATM_0outof1_2(self):
      '''
	   test a pdb file with 1 atom
      0/1 will be selected due to illegal basis_filter
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = 'abc'
      #
      expecting_error = True
      expectd_mask = []
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)

   def test_1ATM_0outof1_3(self):
      '''
	   test a pdb file with 1 atom
      0/1 will be selected due to wrong basis_filter
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = True
      expectd_mask = [0]
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1ATM_1outof1_1(self):
      '''
	   test a pdb file with 1 atom
      1/1 selected by basis_filter ('name[i]=="N"')
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False
      expectd_mask = [1]
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1ATM_1outof1_2(self):
      '''
	   test a pdb file with 1 atom
      1/1 selected by basis_filter ('resid[i]==515')
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False
      expectd_mask = [1]
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1ATM_1outof1_3(self):
      '''
	   test a pdb file with 1 atom
      1/1 selected by basis_filter ('name[i]=="N" and resid[i]==515')
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      basis_filter = 'name[i]=="N" and resid[i]==515'
      #
      expecting_error = False
      expectd_mask = [1]
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_0outof15_1(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      0/15 will be selected due to empty basis_filter
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = ''
      #
      expecting_error = True
      expectd_mask = []
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_0outof15_2(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      0/15 will be selected due to illegal basis_filter
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'abc'
      #
      expecting_error = True
      expectd_mask = []
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_0outof15_3(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      0/15 will be selected due to wrong selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'name[i]=="B" or resid[i]==12'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = True
      expectd_mask = [0]*15
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_8outof15_1(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      8/15 will be selected due to right atom selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O"'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*4+[0]*4+[1]*4+[0]*3
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_8outof15_2(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      8/15 will be selected due to right residue selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'resid[i]==515'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*8+[0]*7
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_2AAD_4outof15_3(self):
      '''
	   test a pdb file with 15 atoms and 2 residue
      3/15 will be selected due to right comprehensive selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'resid[i]==515 and (name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O") and beta[i]>10.0'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*4+[0]*11
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_rna_1(self):
      '''
      test a rna
      0/10632 will be selected due to empty basis_filter
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      basis_filter = ''
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = True
      expectd_mask = [] #[0]*10632
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_rna_2(self):
      '''
      test a rna
      0/10632 will be selected due to wrong basis_filter
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      basis_filter = 'moltype[i]=="protein"'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = True
      expectd_mask = [0]*10632
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_rna_3(self):
      '''
      test a rna
      10632/10632 will be selected due to wrong basis_filter
      This test is hardwired to pass now due to an unfixed bug for the moltype determination between rna and dna
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      #basis_filter = 'moltype[i]=="nucleic"'
      basis_filter = 'moltype[i]=="rna"'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      print error, mask
      #
      expecting_error = False
      expectd_mask = [0]*10632
      #
      self.assertEqual(len(error)>0, expecting_error)
      #self.assertEqual(list(mask),expectd_mask)


   def test_1CRN_1(self):
      '''
	   test a small protein
      0/327 will be selected due to wrong selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      basis_filter = 'chain[i]=="C" and resid[i]==-1 and (name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O")'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = True
      expectd_mask = [0]*(327-0)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1CRN_2(self):
      '''
	   test a small protein
      1/327 will be selected 
      '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      basis_filter = 'chain[i]=="A" and resid[i]==1 and (name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O")'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*4+[0]*(327-4)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1CRN_2(self):
      '''
	   test a small protein
      327/327 will be selected 
      '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      basis_filter = 'chain[i]=="A"'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*327
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files")
   def test_1KP8_1(self):
      '''
	   test a groel
      0/57085 will be selected due to wrong selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      basis_filter = 'chain[i]=="A" and resid[i]==-1 and (name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O")'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = True
      expectd_mask = [0]*(57085-0)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files")
   def test_1KP8_2(self):
      '''
	   test a groel
      0/57085 will be selected due to right atom/residue selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      basis_filter = 'chain[i]=="A" and resid[i]==2 and (name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O")'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*4+[0]*(57085-4)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files")
   def test_1KP8_3(self):
      '''
	   test a groel
      57085/57085 will be selected due to right atom/residue selection
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      basis_filter = 'chain[i]!="Z" and resid[i]<100000'
      #
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      expecting_error = False
      expectd_mask = [1]*(57085-0)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1PSI(self):
      '''
      test a pdb which will not be read successfully by read_pdb
      assertRaises
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      basis_filter = 'name[i]=="N"'
      #
      with self.assertRaises(Exception):
         self.o.get_subset_mask(basis_filter)



   def tearDown(self):
      pass


if __name__ == '__main__': 
   main() 
