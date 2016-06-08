'''
contract

make sure the keys are unique
make sure the right amu list was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol


import os
DataPath = os.path.dirname(os.path.realpath(__file__))+'/../data/sasmol/sasproperties/'


class Test_unit_sasproperties_Atomic_amu(MockerTestCase):

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

