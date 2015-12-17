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

3 residues, mask none due to empty list
3 residues, mask none due to wrong flexible_residues list
3 residues, mask none due to illegal flexible_residues list
3 residues, mask 1st
3 residues, mask 2nd
3 residues, mask 3rd
3 residues, mask 1st and 2nd
3 residues, mask 1st and 3rd
3 residues, mask 3rd and 2nd
3 residues, mask all
3 residues with duplicate resid in 2 chains, mask 1st (both chains will be masked, what about if we only need to mask one?)

small protein (crambin with 46aa), randomly mask [12, 36, 46, 18, 8]
small protein (crambin with 46aa), mask all

large protein complex (groel with 526*14 residues), randomly mask [12, 36, 46, 18, 8] (Skipped as SASSIE_LARGETEST)
large protein complex (groel with 526*14 residues), mask all (Skipped as SASSIE_LARGETEST)

rna molecule, mask randomly 3 residue dihedrals
rna molecule, mask all residues

problemetic pdb (1PSI wih unpaird MODEL/ENDMDL)
"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset
import numpy

import os

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_intg_sassubset_Mask_get_dihedral_subset_mask(MockerTestCase): 
 
   def setUp(self):
      self.o=sasmol.SasMol(0)

   def assert_list_almost_equal(self,a,b,places=5):
        if (len(a)!=len(b)):
           raise TypeError
        else:
           for i in range(len(a)):
              if isinstance(a[i],(int,float,numpy.generic)):
                 if (numpy.isnan(a[i]) and numpy.isnan(b[i])): continue
                 self.assertAlmostEqual(a[i],b[i],places)
              else:
                 self.assert_list_almost_equal(a[i],b[i],places)

   def hard_wired_get_single_dihedral_subset_mask(self, o, nfr, mtype):
      farray = [0]*len(o.name())
      for i in range(len(o.name())):
         name = o.name()[i]
         resid = o.resid()[i]
         if mtype==0:
            if (resid==nfr-1 and name=='C') or (resid == nfr and name in ['N','CA','C']) or (resid==nfr+1 and name=='N'):
               farray[i]=1
         elif mtype==1:
            if (resid==nfr-1 and name=="O3'") or (resid == nfr and name in ["P","O5'","C5'","C4'","C3'","O3'"]) or (resid==nfr+1 and name=='P') or (resid==nfr+1 and name=="O5'"):
               farray[i]=1
      return farray

   def hard_wired_get_all_dihedral_subset_mask(self, o, flexible_residues, mtype):
      farray = []
      for nfr in flexible_residues:
         farray.append(self.hard_wired_get_single_dihedral_subset_mask(o, nfr, mtype))
      return farray
       

   def test_null(self):
      '''
      null test
      read a nonexisting pdb
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      #
      flexible_residues = []
      mtype=0
      with self.assertRaises(AttributeError):
         result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)


   def test_3AAD_mask_none_empty_flexible_residues(self):
      '''
	   test a pdb file with 3 residue
      mask none of the residue dihedral due to empty flexible residue list
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = []
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_none_wrong_flexible_residues(self):
      '''
	   test a pdb file with 3 residue
      mask none of the residue dihedral due to wrong flexible residue list
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [1,2]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_none_illegal_flexible_residues(self):
      '''
	   test a pdb file with 3 residue
      mask none of the residue dihedral due to wrong flexible residue list
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [1,'a']
      mtype=0
      with self.assertRaises(ValueError):
         result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)


   def test_3AAD_mask_1st(self):
      '''
	   test a pdb file with 3 residue
      mask the 1st residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [515]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_2nd(self):
      '''
	   test a pdb file with 3 residue
      mask the 2nd residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [516]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_1st_2nd(self):
      '''
	   test a pdb file with 3 residue
      mask the 1st and 2nd residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [515,516]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_1st_3nd(self):
      '''
	   test a pdb file with 3 residue
      mask the 1st and 3rd residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [515,517]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)



   def test_3AAD_mask_3rd_2nd(self):
      '''
	   test a pdb file with 3 residue
      mask the 3rd and 2nd residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [515,516]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_mask_all(self):
      '''
	   test a pdb file with 3 residue
      mask the all residue dihedrals
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD.pdb')
      #
      flexible_residues = [515,516,517]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_3AAD_2chains_mask_1st(self):
      '''
	   test a pdb file with 3 residues with duplicating resid in 2 chains
      mask the 1st residue dihedral
	   '''
      #
      self.o.read_pdb(PdbDataPath+'3AAD-2chain.pdb')
      #
      flexible_residues = [515]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   def test_1CRN_mask_random_5(self):
      '''
      test a small protein (crambin)
      mask randomly 5 residue dihedrals
      '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      flexible_residues = [12, 36, 46, 18, 8]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)
      
      
   def test_1CRN_mask_all(self):
      '''
      test a small protein (crambin)
      mask the all residue dihedrals
      '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      flexible_residues = range(1,46)
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files")
   def test_1KP8_mask_random_5(self):
      '''
	   test a large protein complex (groel)
      mask randomly 5 residue dihedrals
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      flexible_residues = [12, 36, 46, 18, 8]
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)
      

   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files. It will take 5 min")
   def test_1KP8_mask_all(self):
      '''
	   test a large protein complex (groel)
      mask the all residue dihedrals
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      flexible_residues = range(2,527)
      mtype=0
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)

   def test_rna_mask_random(self):
      '''
      test a rna molecule
      mask randomly 3 residue dihedrals
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      flexible_residues = [1, 7, 13]
      mtype=1
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)

   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing huge files. It will take 5 min")
   def test_rna_mask_all(self):
      '''
      test a rna molecule
      mask all residue dihedrals
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      flexible_residues = range(1,25)
      mtype=1
      result_farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      expected_farray =  self.hard_wired_get_all_dihedral_subset_mask(self.o, flexible_residues, mtype)
      #
      print 'result_mask:\n', result_farray, '\nexpected_mask:\n',expected_farray
      self.assertTrue(isinstance(result_farray, numpy.ndarray))
      self.assert_list_almost_equal(result_farray, expected_farray)


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
      flexible_residues = range(2,527)
      mtype=0
      #
      with self.assertRaises(Exception):
         self.o.get_dihedral_subset_mask(flexible_residues,mtype)


   def tearDown(self):
      pass



if __name__ == '__main__': 
   main() 

