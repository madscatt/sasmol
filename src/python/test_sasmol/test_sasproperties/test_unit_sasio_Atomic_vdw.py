'''
contract

make sure the keys are unique
make sure the right vdw list was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol


import os
DataPath = os.path.dirname(os.path.realpath(__file__))+'/../data/sasmol/sasproperties/'


class Test_unit_sasproperties_Atomic_van_der_Waals_radii(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.vdw  = self.o.van_der_Waals_radii()


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
      vdw_unique = self.unique(self.vdw.keys())
      self.assertEqual(self.vdw.keys(), vdw_unique)



   def test_all(self):
      '''
      make sure the right vdw list was generated
      '''
      #
      datafile = DataPath+'vdw.txt'
      fp = open(datafile,'r')
      ele = {}
      for line in fp.readlines():
         van_der_Waals_radii = line.split()
#         print van_der_Waals_radii
         if van_der_Waals_radii[1] == 'None':
            ele[van_der_Waals_radii[0]] = None
         else:         
            ele[van_der_Waals_radii[0]]=float(van_der_Waals_radii[1])
#         print ele
      self.assertEqual(ele,self.vdw)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

