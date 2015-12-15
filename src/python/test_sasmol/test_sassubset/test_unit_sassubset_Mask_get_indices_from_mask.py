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
contract:

null test, 

test for 1 atom, with none masked
test for 1 atom, with 1 atom masked
test for 1 atom, with none masked
test for 1 atom, with the first atom masked
test for 1 atom, with the second atom masked
test for 1 atom, with both atoms masked
test for 1000 atom, with none masked
test for 1000 atom, with even atomic number masked
test for 1000 atom, with all masked

negative test, with mask length longer than natoms
negative test, with mask length shorter than natoms
negative test, with mask length shorter than natoms

"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset

import numpy

import os

class Test_unit_sassubset_Mask_get_indices_from_mask(MockerTestCase): 

   def mock_up(self, Cls_ptch, mthd, mocker, result=None, mmin=0, mmax=None):
      methodToCall = getattr(Cls_ptch,mthd)
      methodToCall(ARGS)
      mocker.result(result)
      mocker.count(mmin, mmax)

   def mock_up_get_indices_from_mask(self, Cls, mocker, natoms):
      Cls_ptch = mocker.patch(Cls)
      self.mock_up(Cls_ptch, 'natoms', mocker, natoms)
      mocker.replay()

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
 

   def setUp(self):
      self.m = Mocker()
      self.o=sasmol.SasMol(0)


   def test_null(self):
      '''
      null test, 
      '''
      natoms=0
      mask=[]
      #
      expected_indices = []
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1atom_mask_none(self):
      '''
      test for 1 atom, with none masked
      '''
      natoms=1
      mask=[0]
      #
      expected_indices = []
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1atom_mask_one(self):
      '''
      test for 1 atom, with 1 atom masked
      '''
      natoms=1
      mask=[1]
      #
      expected_indices = [0]
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2atom_mask_none(self):
      '''
      test for 1 atom, with none masked
      '''
      natoms=2
      mask=[0,0]
      #
      expected_indices = []
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2atom_mask_first(self):
      '''
      test for 1 atom, with the first atom masked
      '''
      natoms=2
      mask=[1,0]
      #
      expected_indices = [0]
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2atom_mask_second(self):
      '''
      test for 1 atom, with the second atom masked
      '''
      natoms=2
      mask=[0,1]
      #
      expected_indices = [1]
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_2atom_mask_both(self):
      '''
      test for 1 atom, with both atoms masked
      '''
      natoms=2
      mask=[1,1]
      #
      expected_indices = [0,1]
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1000atoms_mask_none(self):
      '''
      test for 1000 atom, with none masked
      '''
      natoms=1000
      mask=[0]*natoms
      #
      expected_indices = []
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1000atoms_mask_even(self):
      '''
      test for 1000 atom, with even atomic number masked
      '''
      natoms=1000
      mask=[1,0]*(natoms/2)
      #
      expected_indices = range(0,natoms,2)
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_1000atoms_mask_all(self):
      '''
      test for 1000 atom, with all masked
      '''
      natoms=1000
      mask=[1]*natoms
      #
      expected_indices = range(natoms)
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      result_indices = self.o.get_indices_from_mask(mask)
      #
      print 'expected_indices:\n', expected_indices, '\nresult_indices:\n',result_indices
      self.assert_list_almost_equal(expected_indices, result_indices)


   def test_negative_1(self):
      '''
      negative test, with mask length longer than natoms
      '''
      natoms=10
      mask=[1]*natoms*2
      #
      expected_indices = range(natoms)
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      with self.assertRaises(Exception):
         result_indices = self.o.get_indices_from_mask(mask)
   

   def test_negative_2(self):
      '''
      negative test, with mask length shorter than natoms
      '''
      natoms=10
      mask=[1]*(natoms/2)
      #
      expected_indices = range(natoms)
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      with self.assertRaises(Exception):
         result_indices = self.o.get_indices_from_mask(mask)

   def test_negative_2(self):
      '''
      negative test, with mask length shorter than natoms
      '''
      natoms=10
      mask=[1]*(natoms/2)
      #
      expected_indices = range(natoms)
      #
      self.mock_up_get_indices_from_mask(self.o, self.m, natoms)
      with self.assertRaises(Exception):
         result_indices = self.o.get_indices_from_mask(mask)

   def tearDown(self):
      self.m.verify()

if __name__ == '__main__': 
   main() 
