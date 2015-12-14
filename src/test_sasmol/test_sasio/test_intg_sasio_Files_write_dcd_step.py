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
import sasmol.dcdio as dcdio

import numpy, os

floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_write_dcd_step(MockerTestCase):

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

   def test_1ATM_frame1(self):
      '''
	   test a dcd with 2 frames based on a 1-atom pdb
      write the 1st frame
	   '''
      #
      pdbFileName = DataPath+'1ATM-1to2.pdb'
      dcdFileName = moduleDataPath+'test-results/1ATM-1to2-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\ncoor\n',result_coor
      expected_coor = numpy.array([[[76.944, 41.799, 41.652]]],floattype)
      sum_expected_coor = 160.395
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_1ATM_frame2(self):
      '''
	   test a dcd with 2 frames based on a 1-atom pdb
      write the 2nd frame
	   '''
      #
      pdbFileName = DataPath+'1ATM-1to2.pdb'
      dcdFileName = moduleDataPath+'test-results/1ATM-1to2-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 1
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      expected_coor = numpy.array([[[73.944, 38.799, 41.652]]],floattype)
      sum_expected_coor = 154.395
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_2AAD_frame1(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
      write the 1st frame
	   '''
      #
      pdbFileName = DataPath+'2AAD-1to3.pdb'
      dcdFileName = moduleDataPath+'test-results/2AAD-1to3-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor\n',result_coor
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]]],floattype)
      sum_expected_coor = 2407.676
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_2AAD_frame2(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
      write the 2nd frame
	   '''
      #
      pdbFileName = DataPath+'2AAD-1to3.pdb'
      dcdFileName = moduleDataPath+'test-results/2AAD-1to3-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 1
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      expected_coor = numpy.array([[[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]]],floattype)
      sum_expected_coor = 140.846
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)




   def test_2AAD_frame3(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
      write the 3rd frame
	   '''
      #
      pdbFileName = DataPath+'2AAD-1to3.pdb'
      dcdFileName = moduleDataPath+'test-results/2AAD-1to3-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 2
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      expected_coor = numpy.array([[[ 73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      sum_expected_coor = 1095.772
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)




   def test_rna_frame1to10_frame1(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 1st frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-5.798, 13.082, 22.068],floattype)
      sum_expected_coor = -43307.390
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_rna_frame1to10_frame6(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 6th frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 5
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.348, 14.130, 20.916],floattype)
      sum_expected_coor = -42996.475
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_rna_frame1to10_frame10(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 10th frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 9
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.392, 14.348, 20.914],floattype)
      sum_expected_coor = -42837.531
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)




   def test_rna_frame1to10_frame1(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 1st frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-5.798, 13.082, 22.068],floattype)
      sum_expected_coor = -43307.390
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_rna_frame1to10_frame6(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 6th frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 5
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.348, 14.130, 20.916],floattype)
      sum_expected_coor = -42996.475
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)



   def test_rna_frame1to10_frame10(self):
      '''
      test a dcd  based on a rna molecule of 10 frames
      write the 10th frame
      '''
      #
      #      #
      pdbFileName = DataPath+'rna-1to10.pdb'
      dcdFileName = moduleDataPath+'test-results/rna-1to10-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 9
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.392, 14.348, 20.914],floattype)
      sum_expected_coor = -42837.531
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)


   def test_1CRN(self):
      '''
      test a dcd based on a small protein (crambin)
      '''
      #
      pdbFileName = DataPath+'1CRN.pdb'
      dcdFileName = moduleDataPath+'test-results/1CRN-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      os.remove(dcdFileName)
      self.assertEqual(len(result_coor[0]),327)
      expected_sample_coor = numpy.array([9.555, 2.856, 3.730],floattype)
      sum_expected_coor = 8509.587
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][100], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8(self):
      '''
      test a dcd based on a large pdb complex (groel)
      '''
      #
      pdbFileName = DataPath+'1KP8.pdb'
      dcdFileName = moduleDataPath+'test-results/1KP8-writedcd_step-test.dcd'
      #
      nnatoms = 1
      frame = 0
      step = 1
      self.o.read_pdb(pdbFileName)
      #fp=dcdio.open_dcd_write(dcdFileName)
      #self.o.write_dcd_header(fp, nnatoms)
      fp = self.o.open_dcd_write(dcdFileName)
      self.o.write_dcd_step(fp, frame, step)
      dcdio.close_dcd_write(fp)
      #
      otest = sasmol.SasMol(0)
      otest.read_dcd(dcdFileName)
      result_coor = otest.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor
      self.assertEqual(len(result_coor[0]),57085)
      expected_sample_coor = numpy.array([104.484, 26.915, 12.538],floattype)
      sum_expected_coor = 6269170.260
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][100], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)
      os.remove(dcdFileName)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

