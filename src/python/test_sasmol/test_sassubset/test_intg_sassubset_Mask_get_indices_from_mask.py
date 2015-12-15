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
Contract of Integration test for Mask.set_descriptor_using_mask(

null test by passing an empty mask

negative test by passing a wrong mask
negative test by passing a mask with wrong dimension

1-atom, mask no atom
1-atom, mask 1 atom

2-aa, mask none
2-aa, mask 2 atoms
2-aa, mask 1 reside
2-aa, mask all atoms

rna, mask no atom
rna, mask 300 atoms
rna, mask all atoms

small protein (crambin), mask no atom
small protein (crambin), mask 46 atoms, set charge
small protein (crambin), mask all atoms

large protein (groel), mask no atom
large protein (groel), mask 100 atoms
large protein (groel), mask all atoms, set moltype (Skipped as SASSIE_LARGETEST)

a bad pdb, assertRaises

"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset

import numpy

import os, copy

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sassubset_Mask_set_descriptor_using_mask(MockerTestCase): 
 

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.o_result=sasmol.SasMol(1)
      self.o_expected=sasmol.SasMol(2)


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


   def hard_wired_set_expected_descriptor(self, mask, d, value):
      for i in range(len(mask)):
         if mask[i]==1:
            d[i]=value
            

   def test_null(self):
      '''
      null test by passing an empty mask
      '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      mask = [0]
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_negative_1(self):
      '''
	   negative test by passing a wrong mask
      '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      mask = 'a'
      #
      with self.assertRaises(Exception):
         result_indices = self.o.get_indices_from_mask(mask)


   def test_negative_2(self):
      '''
	   negative test by passing a mask with wrong dimension
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      mask = [1,0,0]
      #
      with self.assertRaises(Exception):
         result_indices = self.o.get_indices_from_mask(mask)


   def test_1ATM_mask_none(self):
      '''
	   1-atom, mask no atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      mask = [0]
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1ATM_mask_one(self):
      '''
	   1-atom, mask 1 atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      mask = [1]
      #
      expected_indices = [0]
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2AAD_mask_atomNone(self):
      '''
	   2-aa, mask none
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      mask = [0]*15
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2AAD_1frame_mask_2atoms(self):
      '''
	   2-aa, mask 2 atoms
      '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      mask = [1,1]+[0]*13
      #
      expected_indices = [0,1]
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2AAD_1frame_mask_1res(self):
      '''
	   2-aa, mask 1 reside
      '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      mask = [0]*8+[1]*7
      #
      expected_indices = range(8,15)
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2AAD_mask_atoms_all(self):
      '''
	   2-aa, mask all atoms
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      mask = [1]*15
      #
      expected_indices = range(15)
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_rna_1(self):
      '''
      rna, mask no atom
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      mask = [0]*self.o.natoms()
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_rna_2(self):
      '''
      rna, mask 300 atoms
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      mask = [0]*100+[1]*200+[0]*(self.o.natoms()-300)
      #
      expected_indices = range(100,300)
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_rna_3(self):
      '''
      rna, mask all atoms
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      mask = [1]*self.o.natoms()
      #
      expected_indices = range(self.o.natoms())
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1CRN_1(self):
      '''
      small protein (crambin), mask no atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      mask = [0]*self.o.natoms()
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)
      

   def test_1CRN_2(self):
      '''
      small protein (crambin), mask 46 atoms, set charge
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      mask = [0]*20+[1]*30+[0]*10+[1]*(self.o.natoms()-60)
      #
      expected_indices = range(20,50)+range(60,self.o.natoms())
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)



   def test_1CRN_3(self):
      '''
      small protein (crambin), mask all atoms
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      mask = [1]*self.o.natoms()
      #
      expected_indices = range(self.o.natoms())
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)



   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_1(self):
      '''
      large protein (groel), mask no atom
      '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      mask = [0]*self.o.natoms()
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_2(self):
      '''
      large protein (groel), mask 100 atoms
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      mask = [0]*100+[1]*(self.o.natoms()-100)
      #
      expected_indices = range(100,self.o.natoms())
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)



   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      large protein (groel), mask all atoms, set moltype (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      mask = [1]*self.o.natoms()
      #
      expected_indices = range(self.o.natoms())
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)



   def test_1PSI(self):
      '''
      a bad pdb, assertRaises
      '''
      try:
         self.o.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      mask = [0]*self.o.natoms()
      #
      expected_indices = []
      #
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def tearDown(self):
      pass



if __name__ == '__main__': 
   main() 
