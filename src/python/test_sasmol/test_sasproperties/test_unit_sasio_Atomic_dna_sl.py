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


class Test_unit_sasproperties_Atomic_dna_sl(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.dna_sl  = self.o.dna_sl()


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
      dna_sl_unique = self.unique(self.dna_sl.keys())
      self.assertEqual(self.dna_sl.keys(), dna_sl_unique)


   def test_all(self):
      '''
      make sure the right element_sl list was generated
      '''
      #
      datafile = DataPath+'dna_sl.txt'
      fp = open(datafile,'r')
      ele = {}
      for line in fp.readlines():
         residue_scattering = line.split()
#         print residue_scattering       
         base = residue_scattering[0]
#         print base
         value = eval(line[4:])
#         print value
         ele[base]=value
#         print ele
      self.assertEqual(ele,self.dna_sl)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

