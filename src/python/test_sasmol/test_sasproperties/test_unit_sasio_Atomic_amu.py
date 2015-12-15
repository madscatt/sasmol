
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


contract

make sure the keys are unique
make sure the right amu list was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasproperties')+os.path.sep

class Test_unit_sasproperties_Atomic_charmm_names(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.standard_atomic_weight  = self.o.amu()


   def unique(self,seq):
      '''
      remove the duplicate elements
      '''
      seen = set()
      seen_add = seen.add
      return [ x for x in seq if x not in seen and not seen_add(x)]


   def test_uniqe(self):
      '''
	   make sure the keys are unique
      '''
      #
      standard_atomic_weight_unique = self.unique(self.standard_atomic_weight.keys())
      self.assertEqual(self.standard_atomic_weight.keys(), standard_atomic_weight_unique)



   def test_all(self):
      '''
      make sure the right amu list was generated
      '''
      #
      datafile = DataPath+'standard_atomic_weight.txt'
      fp = open(datafile,'r')
      ele = {}
      for line in fp.readlines():
         amu = line.split()
         ele[amu[0]]=float(amu[1])
      self.assertEqual(ele,self.standard_atomic_weight)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

