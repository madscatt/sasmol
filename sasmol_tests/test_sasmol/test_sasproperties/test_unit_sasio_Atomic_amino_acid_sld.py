'''
contract

make sure the keys are unique
make sure the right amu list was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

from sassie.sasmol import sasmol


import os
DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/sasmol/sasproperties/'


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

