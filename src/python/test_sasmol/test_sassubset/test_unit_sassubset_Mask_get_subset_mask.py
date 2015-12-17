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

'''
unit test for get_subset_mask:

contract

1. null list
2. 0/1 masked due to empty basis_filter
3. 0/1 masked due to wrong basis_filter
4. 0/1 masked due to wrong selection
5. 1/1 masked due to right atom selection
6. 1/1 masked due to right residue selection
7. 1/1 masked due to right residue/atom selection
8. 3 atoms...
9. 500 atoms...
10. 5000000 atoms...
'''

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset
import numpy

import os

class Test_sassubset_Mask_get_subset_mask(MockerTestCase): 

   """
   #NOT SURE WHY DOESNT WORK
   def mock_it_up(self, mocker, function, result=None, mmin=0, mmax=None):
      function = mocker.mock()
      function(ARGS)
      mocker.result(result)
      mocker.count(mmin, mmax)
   """

   def set_parameters(self, index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms):
      self.o.setIndex(index)
      self.o.setName(name)
      self.o.setLoc(loc)
      self.o.setResname(resname)
      self.o.setChain(chain)
      self.o.setResid(resid)
      self.o.setRescode(rescode)
      self.o.setOccupancy(occupancy)
      self.o.setBeta(beta)
      self.o.setSegname(segname)
      self.o.setElement(element)
      self.o.setCharge(charge)
      self.o.setMoltype(moltype)
      self.o.setNatoms(natoms)
 

   def setUp(self):
      self.o=sasmol.SasMol(0)


   def test_null(self):
      '''
      test for a null atom list
      '''
      #
      index=[]
      name=[]
      loc=[]
      resname=[]
      chain=[]
      resid=[]
      rescode=[]
      occupancy=[]
      beta=[]
      segname=[]
      element=[]
      charge=[]
      moltype=[]
      natoms=0
      #
      basis_filter = ''
      #
      expecting_error = True
      expectd_mask = []
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_single_unmasked_1(self):
      '''
      test for for single atom which will not be masked
      '''
      natoms=1
      index=range(1,natoms+1)
      name=['CA']*natoms
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=range(1,natoms+1)
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C']*natoms
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'resid[i]==3'
      #
      expecting_error = True
      expectd_mask = [0]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_single_umasked_2(self):
      '''
      test for for single atom which will be masked
      '''
      natoms=1
      index=range(1,natoms+1)
      name=['CA']*natoms
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=range(1,natoms+1)
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C']*natoms
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA" and resid[i]=2'
      #
      expecting_error = True
      expectd_mask = []
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_single_unmasked_3(self):
      '''
      test for for single atom which will not be masked due to filter error
      '''
      natoms=1
      index=range(1,natoms+1)
      name=['CA']*natoms
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=range(1,natoms+1)
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C']*natoms
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'resid[js]==3'
      #
      expecting_error = True
      expectd_mask = []
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_single_masked_1(self):
      '''
      test for for single atom which will be masked
      '''
      natoms=1
      index=range(1,natoms+1)
      name=['CA']*natoms
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=range(1,natoms+1)
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C']*natoms
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #      
      basis_filter = 'name[i]=="CA"'
      #
      expecting_error = False
      expectd_mask = [1]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_single_masked_2(self):
      '''
      test for for single atom which will be masked
      '''
      natoms=1
      index=range(1,natoms+1)
      name=['CA']*natoms
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=range(1,natoms+1)
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C']*natoms
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA" and resid[i]==1'
      #
      expecting_error = False
      expectd_mask = [1]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_0out3_masked(self):
      '''
      test for for three atoms where none will be masked
      '''
      natoms=3
      index=range(1,natoms+1)
      name=['CA','N','C']
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C']
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = True
      expectd_mask = [0,0,0]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1out3_masked_1(self):
      '''
      test for for three atoms where one will be masked
      '''
      natoms=3
      index=range(1,natoms+1)
      name=['CA','N','C']
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C']
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA"'
      expecting_error = False
      expectd_mask = [1,0,0]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_1out3_masked_2(self):
      '''
      test for for three atoms where one will be masked
      '''
      natoms=3
      index=range(1,natoms+1)
      name=['CA','N','C']
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C']
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="N" and resid[i]==1'
      expecting_error = False
      expectd_mask = [0,1,0]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_3out3_masked_1(self):
      '''
      test for for three atoms where all will be masked
      '''
      natoms=3
      index=range(1,natoms+1)
      name=['CA','N','C']
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C']
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'resid[i]==1'
      #
      expecting_error = False
      expectd_mask = [1,1,1]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_3out3_masked(self):
      '''
      test for for three atoms where none will be masked
      '''
      natoms=3
      index=range(1,natoms+1)
      name=['CA','N','C']
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C']
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'resid[i]==1'
      #
      expecting_error = False
      expectd_mask = [1,1,1]
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expectd_mask)


   def test_100out500_masked(self):
      '''
      test for 500 atoms where 100 will be masked
      '''
      natoms=500
      aasize=5
      index=range(1,natoms+1)
      name=['N','CA','C','O','CB']*(natoms/aasize)
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C','O','C']*(natoms/aasize)
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA"'
      #
      expecting_error = False      
      expected_mask = [0,1,0,0,0]*(natoms/aasize)
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expected_mask)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_10000out50000_masked(self):
      '''
      test for 50000 atoms where 10000 will be masked
      '''
      natoms=50000
      aasize=5
      index=range(1,natoms+1)
      name=['N','CA','C','O','CB']*(natoms/aasize)
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C','O','C']*(natoms/aasize)
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA"'
      #
      expecting_error = False      
      expected_mask = [0,1,0,0,0]*(natoms/aasize)
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expected_mask)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
   def test_1000000out5000000_masked(self):
      '''
      test for 5000000 atoms where 1000000 will be masked
      '''
      natoms=5000000
      aasize=5
      index=range(1,natoms+1)
      name=['N','CA','C','O','CB']*(natoms/aasize)
      loc=[' ']*natoms
      resname=['ALA']*natoms
      chain=['A']*natoms
      resid=[1]*natoms
      rescode=[' ']*natoms
      occupancy=[1.0]*natoms
      beta=[10.00]*natoms
      segname=[' ']*natoms
      element=['C','N','C','O','C']*(natoms/aasize)
      charge=[1.0]*natoms
      moltype=['protein']*natoms
      #
      basis_filter = 'name[i]=="CA"'
      #
      expecting_error = False      
      expected_mask = [0,1,0,0,0]*(natoms/aasize)
      #
      self.set_parameters(index, name, loc, resname, chain, resid, rescode, occupancy, beta, segname, element, charge, moltype, natoms)
      error, mask = self.o.get_subset_mask(basis_filter)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assertEqual(list(mask),expected_mask)






   def tearDown(self):
      pass

if __name__ == '__main__': 
   main() 

