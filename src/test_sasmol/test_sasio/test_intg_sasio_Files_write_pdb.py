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

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop
import sasmol.sascalc as sascalc

import numpy, os, copy

import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep


class Test_intg_sasio_Files_write_pdb(MockerTestCase):

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
	   test a 1-atom pdb
	   '''
      #
      frame = 0
      self.o.read_pdb(DataPath+'1ATM.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/1ATM-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1ATM-writepdb-test.pdb')


   def test_1ATM_frame1to2_second_frame(self):
      '''
	   test the second frame a pdb file with 1 atom and 2 frames
	   '''
      #
      frame = 1
      self.o.read_pdb(DataPath+'1ATM-1to2.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/1ATM-1to2-2-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1ATM-1to2-2-writepdb-test.pdb')


   def test_2AAD(self):
      '''
	   test a 2-aa pdb
	   '''
      #
      frame = 0
      self.o.read_pdb(DataPath+'2AAD.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/2AAD-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/2AAD-writepdb-test.pdb')


   def test_2AAD_frames1to3_second_frame(self):
      '''
	   test the second frame a pdb file with 2 amino acids and 3 frames
	   '''
      #
      frame = 2
      self.o.read_pdb(DataPath+'2AAD-1to3.pdb')
      result = self.o.write_pdb(moduleDataPath+'test-results/2AAD-1to3-2-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/2AAD-1to3-2-writepdb-test.pdb')

   def test_rna_frame1to10_frame_3(self):
      '''
	   test the 10th frame of a pdb file of rna with 10 frames
      '''
      #
      frame = 3
      self.o.read_pdb(DataPath+"rna-1to10.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/rna-1to10-3-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/rna-1to10-3-writepdb-test.pdb')


   def test_1CRN(self):
      '''
	   test a small protein (crambin)
      '''
      #
      frame = 0
      self.o.read_pdb(DataPath+"1CRN.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/1CRN-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1CRN-writepdb-test.pdb')


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")   
   def test_1KP8(self):
      '''
	   test a large protein complex (groel)
      '''
      #
      frame = 0
      self.o.read_pdb(DataPath+"1KP8.pdb")
      result = self.o.write_pdb(moduleDataPath+'test-results/1KP8-writepdb-test.pdb',frame,'w')
      self.assertEqual(result,1)
      os.remove(moduleDataPath+'test-results/1KP8-writepdb-test.pdb')



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

