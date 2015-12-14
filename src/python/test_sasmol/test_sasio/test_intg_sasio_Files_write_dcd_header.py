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

from sasmol.test_sasmol.util import env, util

from unittest import main, skipIf 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol
import sasmol.dcdio as dcdio

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_write_dcd_header(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)

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

   def test_1ATM(self):
      '''
	   test a dcd with 2 frames based on a 1-atom pdb
	   '''
      #
      pdbFileName = DataPath+'1ATM.pdb'
      self.o.read_pdb(pdbFileName)
      nset = 2
      #
      dcdFileName = moduleDataPath+'test-results/1ATM-writedcd_header-test.dcd'
      fp = dcdio.open_dcd_write(dcdFileName)
      self.o.write_dcd_header(fp, nset)
      dcdio.close_dcd_write(fp)
      fsize = os.path.getsize(dcdFileName)
      print fsize
      self.assertTrue(fsize>0)
      os.remove(dcdFileName)


   def test_2AAD(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
	   '''
      #
      pdbFileName = DataPath+'2AAD.pdb'
      self.o.read_pdb(pdbFileName)
      nset = 3
      #
      dcdFileName = moduleDataPath+'test-results/2AAD-writedcd_header-test.dcd'
      fp = dcdio.open_dcd_write(dcdFileName)
      self.o.write_dcd_header(fp, nset)
      dcdio.close_dcd_write(fp)
      fsize = os.path.getsize(dcdFileName)
      print fsize
      self.assertTrue(fsize>0)
      os.remove(dcdFileName)


   def test_1CRN(self):
      '''
	   test a dcd from a small protein (crambin)
	   '''
      #
      pdbFileName = DataPath+'1CRN.pdb'
      self.o.read_pdb(pdbFileName)
      nset = 3
      #
      dcdFileName = moduleDataPath+'test-results/1CRN-writedcd_header-test.dcd'
      fp = dcdio.open_dcd_write(dcdFileName)
      self.o.write_dcd_header(fp, nset)
      dcdio.close_dcd_write(fp)
      fsize = os.path.getsize(dcdFileName)
      print fsize
      self.assertTrue(fsize>0)
      os.remove(dcdFileName)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8(self):
      '''
	   test a dcd from a large protein complex (groel)
	   '''
      #
      pdbFileName = DataPath+'1KP8.pdb'
      self.o.read_pdb(pdbFileName)
      nset = 3
      #
      dcdFileName = moduleDataPath+'test-results/1KP8-writedcd_header-test.dcd'
      fp = dcdio.open_dcd_write(dcdFileName)
      self.o.write_dcd_header(fp, nset)
      dcdio.close_dcd_write(fp)
      fsize = os.path.getsize(dcdFileName)
      print fsize
      self.assertTrue(fsize>0)
      os.remove(dcdFileName)


   def test_rna_1to10frames(self):
      '''
	   test a dcd of 10 frames based on a rna molecule
	   '''
      #
      pdbFileName = DataPath+'rna.pdb'
      self.o.read_pdb(pdbFileName)
      nset = 10
      #
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_header-test.dcd'      
      fp = dcdio.open_dcd_write(dcdFileName)
      self.o.write_dcd_header(fp, nset)
      dcdio.close_dcd_write(fp)
      fsize = os.path.getsize(dcdFileName)
      print fsize
      self.assertTrue(fsize>0)
      os.remove(dcdFileName)




   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

