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
Contract of Integration test for Mask.set_descriptor_using_mask

null test by passing a empty descriptor

negative test by passing a wrong mask
negative test by setting natoms
negative test by setting coor

1-atom, mask no atom, set name
1-atom, mask 1 atom, set name
1-atom, mask 1 atom, set element

2-aa, mask no atom, set index
2-aa, mask 2 atoms, set index
2-aa, mask all atoms, set loc

rna, mask no atom, set resid
rna, mask 248 atoms, set occupancy
rna, mask all atoms, set segname

small protein (crambin), mask no atom, set beta
small protein (crambin), mask 46 atoms, set charge
small protein (crambin), mask all atoms, set moltype

large protein (groel), mask no atom, set atom (Skipped as SASSIE_LARGETEST)
large protein (groel), mask 7350 atoms, set index (Skipped as SASSIE_LARGETEST)
large protein (groel), mask all atoms, set moltype (Skipped as SASSIE_LARGETEST)

a bad pdb, assertRaises (some problem)

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
      null test by passing a empty descriptor
      '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      mask = [1]
      descriptor = ''
      value = 100
      expected_descriptor = value
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)


   def test_negative_0(self):
      '''
	   negative test by passing a wrong mask
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      expecting_error = False 
      mask = 1
      descriptor = 'abc'
      value = 100
      expected_descriptor = value
      #
      with self.assertRaises(Exception):
         error = self.o.set_descriptor_using_mask(mask, descriptor, value)


   def test_negative_1(self):
      '''
	   negative test by setting natoms
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      expecting_error = False 
      mask = [1]
      descriptor = self.o.natoms()
      value = 100
      expected_descriptor = value
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)


   def test_negative_2(self):
      '''
	   negative test by setting coor
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      #
      expecting_error = False 
      descriptor = self.o.coor()
      value = [1.0,2.0,3.0]
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)


   def test_1ATM_mask_none(self):
      '''
	   1-atom, mask no atom, set name
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      expecting_error = False 
      mask = [0]
      descriptor = self.o.name()
      value = 'C'
      expected_descriptor = ['N']
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_1ATM_mask_one_1(self):
      '''
	   1-atom, mask 1 atom, set name
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      expecting_error = False 
      mask = [1]
      descriptor = self.o.name()
      value = 'C'
      expected_descriptor = [value]
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_1ATM_mask_one_2(self):
      '''
	   1-atom, mask 1 atom, set element
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      expecting_error = False 
      mask = [1]
      descriptor = self.o.element()
      value = 'O'
      expected_descriptor = [value]
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_2AAD_mask_atomNone(self):
      '''
	   2-aa, mask no atom, set index
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]=="B"'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.index()
      value = 10
      expected_descriptor = copy.deepcopy(descriptor)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_descriptor, descriptor)


   def test_2AAD_1frame_mask_2atoms(self):
      '''
	   2-aa, mask 2 atoms, set index
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]=="N"'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.index()
      value = 10
      expected_descriptor = [10,  2,  3,  4,  5,  6,  7,  8, 10, 10, 11, 12, 13, 14, 15]
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_descriptor, descriptor)


   def test_2AAD_mask_atoms_all(self):
      '''
	   2-aa, mask all atoms, set loc
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]!="B"'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.index()
      value = 10
      expected_descriptor = [value]*self.o.natoms()
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_descriptor, descriptor)


   def test_rna_1(self):
      '''
      rna, mask no atom, set resid
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]==515'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.resid()
      value = 10
      expected_descriptor = copy.deepcopy(descriptor)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_descriptor, descriptor)


   def test_rna_2(self):
      '''
      rna, mask 248 atoms, set occupancy
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]==20'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.occupancy()
      value = '0.5'
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_rna_3(self):
      '''
      rna, mask all atoms, set segname
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]>0'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.segname()
      value = 'X'
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_1CRN_1(self):
      '''
      small protein (crambin), mask no atom, set beta
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]>100000'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.beta()
      value = 10.0
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)
      

   def test_1CRN_2(self):
      '''
      small protein (crambin), mask 46 atoms, set charge
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]=="N"'    
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.charge()
      value = 1.0
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)



   def test_1CRN_3(self):
      '''
      small protein (crambin), mask all atoms, set moltype
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]<100000'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.name()
      value = 'CA'
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_1(self):
      '''
      large protein (groel), mask no atom, set atom (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]=="B"'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.name()
      value = 'CA'
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_2(self):
      '''
      large protein (groel), mask 7350 atoms, set index (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      expecting_error = False
      basis_filter = 'name[i]=="N"'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.index()
      value = 10
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_list_almost_equal(expected_descriptor, descriptor)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      large protein (groel), mask all atoms, set moltype (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      expecting_error = False
      basis_filter = 'resid[i]<100000'      
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      #
      descriptor = self.o.moltype()
      value = "protein"
      expected_descriptor = copy.deepcopy(descriptor)
      self.hard_wired_set_expected_descriptor(mask, expected_descriptor, value)
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)
      #
      print '\nexpected_descriptor:\n',expected_descriptor,'\ndescriptor\n',descriptor
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(expected_descriptor, descriptor)


   def test_1PSI(self):
      '''
      a bad pdb, assertRaises
      '''
      try:
         self.o.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      expecting_error = False 
      mask = [1,0]
      descriptor = ''
      value = 100
      expected_descriptor = value
      #
      error = self.o.set_descriptor_using_mask(mask, descriptor, value)


   def tearDown(self):
      pass


if __name__ == '__main__': 
   main() 
