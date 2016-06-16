'''
contract

make sure the keys are unique
make sure the right element list was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol


import os
DataPath = os.path.dirname(os.path.realpath(__file__))+'/../data/sasmol/sasproperties/'


class Test_unit_sasproperties_Atomic_element_sl(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.element_sl  = self.o.element_sl()


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
      element_sl_unique = self.unique(self.element_sl.keys())
      self.assertEqual(self.element_sl.keys(), element_sl_unique)


   def test_all(self):
      '''
      make sure the right element_sl list was generated
      '''
      #
      datafile = DataPath+'element_sl.txt'
      fp = open(datafile,'r')
      ele = {}
      for line in fp.readlines():
         atomic_scattering = line.split()
#         print atomic_scattering       
         aa = atomic_scattering[0]
#         print aa
         value = eval(line[4:])
#         print value
         ele[aa]=value
#         print ele
      self.assertEqual(ele,self.element_sl)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

