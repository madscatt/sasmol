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

class Test_unit_sasproperties_Atomic_amino_acid_sld(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.amino_acid_sld  = self.o.amino_acid_sld()


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
      amino_acid_sld_unique = self.unique(self.amino_acid_sld.keys())
      self.assertEqual(self.amino_acid_sld.keys(), amino_acid_sld_unique)


   def test_all(self):
      '''
      make sure the right amino_acid_sld list was generated
      '''
      #
      datafile = DataPath+'amino_acid_sld.txt'
      fp = open(datafile,'r')
      ele = {}
      for line in fp.readlines():
         aa = line[0:3]
         value = eval(line[4:])
         ele[aa]=value
      self.assertEqual(ele,self.amino_acid_sld)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

