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
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import numpy
import os

floattype=os.environ['SASSIE_FLOATTYPE']

pdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
dcdDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_write_dcd(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.prcsn = 3

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
      pdbFile = pdbDataPath+'1ATM-1to2.pdb'
      dcdFile = moduleDataPath+'test-results/1ATM-1to2-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      expected_coor = numpy.array([[[76.944, 41.799, 41.652]],[[73.944, 38.799, 41.652]]],floattype)
      sum_expected_coor = 314.790
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)


   def test_2AAD(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
	   '''
      #
      pdbFile = pdbDataPath+'2AAD-1to3.pdb'
      dcdFile = moduleDataPath+'test-results/2AAD-1to3-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]],\
      [[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]],\
      [[  73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)      
      sum_expected_coor = 3644.294
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)


   def test_rna(self):
      '''
	   test a dcd  with 1 frame based on a rna molecule
	   '''
      #
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = moduleDataPath+'test-results/rna-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      self.assertEqual(len(result_coor),1)
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor = 59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631],self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   def test_rna_frame1to10(self):
      '''
      test a dcd with 10 frames based on a rna molecule
      '''
      #
      pdbFile = pdbDataPath+'rna-1to10.pdb'
      dcdFile = moduleDataPath+'test-results/rna-1to10-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      self.assertEqual(len(result_coor),10)
      self.assertEqual(len(result_coor[9]),10632)
      expected_sample_coor = numpy.array([-6.392, 14.348, 20.914],floattype)
      sum_expected_coor = -430804.378
      self.assert_list_almost_equal(expected_sample_coor,result_coor[9][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   def test_1CRN(self):
      '''
      test a dcd from a small protein of crambin
      '''
      #
      pdbFile = pdbDataPath+'1CRN.pdb'
      dcdFile = moduleDataPath+'test-results/1CRN-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      self.assertEqual(len(result_coor[0]),327)
      expected_sample_coor = numpy.array([9.555, 2.856, 3.730],floattype)
      sum_expected_coor = 8509.587
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][100], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")   
   def test_1KP8(self):
      '''
      test a dcd from a protein complex of greol
      '''
      #
      pdbFile = pdbDataPath+'1KP8.pdb'
      dcdFile = moduleDataPath+'test-results/1KP8-writedcd.dcd'
      self.o.read_pdb(pdbFile)
      self.o.write_dcd(dcdFile)
      o1 = sasmol.SasMol(0)
      o1.read_dcd(dcdFile)
      result_coor = o1.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      os.remove(dcdFile)
      #
      self.assertEqual(len(result_coor[0]),57085)
      expected_sample_coor = numpy.array([104.484, 26.915, 12.538],floattype)
      sum_expected_coor = 6269170.260
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][100], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)




   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

