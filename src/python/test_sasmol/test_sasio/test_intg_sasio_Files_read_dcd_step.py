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

from sasmol.test_sasmol.util import env, util, generate_huge_dcd_onthefly

from unittest import main,skipIf
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import numpy
import os

floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep
pdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_read_dcd_step(MockerTestCase):

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
	   test a dcd with 2 frames based on a 1-atom pdb for the first frame
	   '''
      #
      dcdFile = DataPath+'1ATM.dcd'
      frame = 0
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      expected_coor = numpy.array([[[76.944, 41.799, 41.652]]],floattype)
      sum_expected_coor = sum(sum(sum(expected_coor)))
      print '\nexpected_coor \n',expected_coor,'\nsum of expected_coor\n',sum_expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)

   def test_1ATM_frame2(self):
      '''
	   test a dcd with 2 frames based on a 1-atom pdb for the second frame
	   '''
      #
      dcdFile = DataPath+'1ATM.dcd'
      frame = 1
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      expected_coor = numpy.array([[[73.944, 38.799, 41.652]]],floattype)
      sum_expected_coor = sum(sum(sum(expected_coor)))
      print '\nexpected_coor \n',expected_coor,'\nsum of expected_coor\n',sum_expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)

   def test_2AAD_frame1(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb for the first frame
	   '''
      #
      dcdFile = DataPath+'2AAD.dcd'
      frame = 0
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]]],floattype)
      sum_expected_coor = 2407.676
      print '\nexpected_coor \n',expected_coor,'\nsum of expected_coor\n',sum_expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)


   def test_2AAD_frame2(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb for the second frame
	   '''
      #
      dcdFile = DataPath+'2AAD.dcd'
      frame = 1
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      expected_coor = numpy.array([[[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]]],floattype)
      sum_expected_coor = 140.846
      print '\nexpected_coor \n',expected_coor,'\nsum of expected_coor\n',sum_expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)


   def test_2AAD_frame3(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb for the third frame
	   '''
      #
      dcdFile = DataPath+'2AAD.dcd'
      frame = 2
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      expected_coor = numpy.array([[[ 73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      sum_expected_coor = 1095.772
      print '\nexpected_coor \n',expected_coor,'\nsum of expected_coor\n',sum_expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor, self.prcsn)
      self.assertAlmostEqual(sum_expected_coor, sum_result_coor, self.prcsn)



   def test_rna_frame1to10_frame1(self):
      '''
      test a dcd with 10 frames based on a rna molecule for the 1st frame
      '''
      #
      dcdFile = DataPath+'rna-1to10.dcd'
      frame = 0
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-5.798, 13.082, 22.068],floattype)
      sum_expected_coor = -43307.390
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   def test_rna_frame1to10_frame6(self):
      '''
      test a dcd with 10 frames based on a rna molecule for the 6th frame
      '''
      #
      dcdFile = DataPath+'rna-1to10.dcd'
      frame = 5
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.348, 14.130, 20.916],floattype)
      sum_expected_coor = -42996.475
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)



   def test_rna_frame1to10_frame10(self):
      '''
      test a dcd with 10 frames based on a rna molecule for the 10th frame
      '''
      #
      dcdFile = DataPath+'rna-1to10.dcd'
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-6.392, 14.348, 20.914],floattype)
      sum_expected_coor = -42837.531
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_1point0gb(self):
      '''
      test a dcd 1.0gb based on a rna molecule for the 10th frame
      '''
      #
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-1.0g.dcd'
      frames=7813
      generate_huge_dcd_onthefly.generate(fin=pdbFile,fout=dcdFile,frames=frames)
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor = 59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)



   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_1point2gb(self):
      '''
      test a dcd 1.2gb based on a rna molecule
      '''
      #
      dcdFile = '/tmp/rna-1.2g.dcd'
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor =59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_2point0gb(self):
      '''
      test a dcd 2.0gb based on a rna molecule
      '''
      #
      dcdFile = '/tmp/rna-2.0g.dcd'
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor = 59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_3point2gb(self):
      '''
      test a dcd 3.2gb based on a rna molecule
      '''
      #
      dcdFile = '/tmp/rna-3.2g.dcd'
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor =  59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_6point4gb(self):
      '''
      test a dcd 6.4gb based on a rna molecule
      '''
      #
      dcdFile = '/tmp/rna-6.4g.dcd'
      frame = 9
      fp = self.o.open_dcd_read(dcdFile)
      self.o.setCoor(numpy.zeros((1,fp[1],3),floattype))
      for i in range(frame+1):
         self.o.read_dcd_step(fp,frame)
      result_coor = self.o.coor()
      sum_result_coor = sum(sum(sum(result_coor)))
      print '\nresult_coor \n',result_coor,'\nsum of result_coor\n',sum_result_coor
      #
      self.assertEqual(len(result_coor[0]),10632)
      expected_sample_coor = numpy.array([-35.960, 48.536, 72.994],floattype)
      sum_expected_coor = 59505.827
      self.assert_list_almost_equal(expected_sample_coor,result_coor[0][10631], self.prcsn)
      self.assertAlmostEqual(sum_result_coor, sum_expected_coor, self.prcsn)



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

