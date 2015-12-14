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
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol

import numpy, os, copy

floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_read_dcd(MockerTestCase):

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

   def test_1ATM_one_frame(self):
      '''
	   test a pdb file with 1 atom and 1 frame
	   '''
      #
      self.o.read_pdb(DataPath+'1ATM.pdb')
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor
      #
      expected_coor = numpy.array([[[73.944, 41.799, 41.652]]],floattype)
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)


   def test_1ATM_two_frames(self):
      '''
	   test a pdb file with 1 atom and 2 frames
	   '''
      #
      self.o.read_pdb(DataPath+'1ATM-1to2.pdb')
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor
      #
      expected_coor = numpy.array([[[76.944, 41.799, 41.652]],[[73.944, 38.799, 41.652]]],floattype)
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)


   def test_2AAD_one_frame(self):
      '''
	   test a pdb file with 2 amino acids and 1 frame
	   '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]]],floattype)
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)


   def test_2AAD_three_frames_separatedby_END(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by END)
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-END.pdb')
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]],\
      [[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]],\
      [[  73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)
 

   def test_2AAD_three_frames_separatedby_END_all_properties(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by END) for all properties
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-END.pdb')
      natoms = self.o.natoms()
      self.assertEqual(natoms,15)
      self.assertEqual(self.o.moltype(),['protein']*natoms)
      self.assertEqual(self.o.number_of_frames(),3)
      self.assertEqual(self.o.atom(),['ATOM']*natoms)
      self.assertEqual(self.o.name(),['N','CA','C','O','CB','CG1','CG2','CD1','N','CA','C','O','CB','OG1','CG2'])
      self.assert_list_almost_equal(self.o.index(),range(1,natoms+1))
      self.assertEqual(self.o.loc(),[' ']*natoms)
      self.assertEqual(self.o.resname(),['ILE']*8+['THR']*7)
      self.assertEqual(self.o.chain(),['N']*natoms)
      self.assert_list_almost_equal(self.o.resid(),[515]*8+[516]*7)
      self.assertEqual(self.o.rescode(),[' ']*natoms)
      self.assertEqual(self.o.occupancy(),['1.00']*natoms)
      self.assertEqual(self.o.beta(),['36.37', '36.23', '36.32', '36.04', '36.69', '38.12', '34.42', '39.85', '35.01', '35.51', '38.09', '36.94', '36.74', '37.19', '34.44'])
      self.assertEqual(self.o.segname(),['N']*natoms)
      self.assertEqual(self.o.element(),['N', 'C', 'C', 'O', 'C', 'C', 'C', 'C', 'N', 'C', 'C', 'O', 'C', 'O', 'C'])
      self.assertEqual(self.o.charge(),['  ']*natoms)


   def test_2AAD_three_frames_separatedby_END_wrong_number_atoms(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-END_wrong_number_atoms.pdb')



   def test_2AAD_three_frames_separatedby_MODEL(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL.pdb')
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor
      #
      expected_coor = numpy.array([[[  73.944,   41.799,   41.652], [  74.229,   42.563,   40.456], [  75.667,   43.093,   40.463], [  76.264,   43.279,   39.401], [  73.210,   43.734,   40.336], [  71.856,   43.168,   39.926], [  73.677,   44.782,   39.354], [  70.721,   44.177,   39.946], [  76.231,   43.330,   41.647], [  77.592,   43.852,   41.730], [  78.617,   42.820,   42.184], [  79.712,   43.169,   42.656], [  77.671,   45.097,   42.648], [  77.054,   44.816,   43.910], [  76.970,   46.273,   42.000]],\
      [[ -73.944,   41.799,   41.652], [ -74.229,   42.563,   40.456], [ -75.667,   43.093,   40.463], [ -76.264,   43.279,   39.401], [ -73.210,   43.734,   40.336], [ -71.856,   43.168,   39.926], [ -73.677,   44.782,   39.354], [ -70.721,   44.177,   39.946], [ -76.231,   43.330,   41.647], [ -77.592,   43.852,   41.730], [ -78.617,   42.820,   42.184], [ -79.712,   43.169,   42.656], [ -77.671,   45.097,   42.648], [ -77.054,   44.816,   43.910], [ -76.970,   46.273,   42.000]],\
      [[  73.944,  -41.799,   41.652], [  74.229,  -42.563,   40.456], [  75.667,  -43.093,   40.463], [  76.264,  -43.279,   39.401], [  73.210,  -43.734,   40.336], [  71.856,  -43.168,   39.926], [  73.677,  -44.782,   39.354], [  70.721,  -44.177,   39.946], [  76.231,  -43.330,   41.647], [  77.592,  -43.852,   41.730], [  78.617,  -42.820,   42.184], [  79.712,  -43.169,   42.656], [  77.671,  -45.097,   42.648], [  77.054,  -44.816,   43.910], [  76.970,  -46.273,   42.000]]],floattype)
      print '\nexpected_coor \n',expected_coor
      #
      self.assert_list_almost_equal(expected_coor, result_coor,3)
 

   def test_2AAD_three_frames_separatedby_MODEL_wrong_number_atoms(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):      
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_wrong_number_atoms.pdb')



   def test_2AAD_three_frames_separatedby_MODEL_wrongnumber_mix_END(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_wrongnumber_mix_END.pdb')



   def test_2AAD_three_frames_separatedby_MODEL_mix_END_noterminating(self):
      '''
	   test a pdb file with 2 amino acids and 3 frames (separated by MODEL)
	   '''
      #
      with self.assertRaises(Exception):      
         self.o.read_pdb(moduleDataPath+'2AAD-1to3-MODEL_mix_END_noterminating.pdb')


   
   def test_rna_frame1to10_frame_3(self):
      '''
	   test a pdb file of rna with 10 frames
      '''
      #
      self.o.read_pdb(DataPath+"rna-1to10.pdb")
      result_coor = self.o.coor()
      print '\nresult_coor \n',result_coor[2][10627]
      #
      self.assertEqual(len(result_coor),10)
      self.assertEqual(len(result_coor[3]),10632)
      expected_coor_sample = numpy.array([-5.564, 20.324, 26.185],floattype) #atom 10627 of frame 3
      self.assert_list_almost_equal(result_coor[2][10627],expected_coor_sample,3)


   def test_1PSI(self):
      '''
	   test a pdb file without ENDMDL
	   '''
      #
      with self.assertRaises(Exception):
         self.o.read_pdb(DataPath+'1PSI.pdb')

   def test_blanklines(self):
      '''
	   test a pdb file ending with blank lines
	   '''
      #
      self.o.read_pdb(DataPath+'dimcd_fixed_atoms.pdb')
      expected_coor_sample = numpy.array([65.124,  35.624,  50.733],floattype)
      result_coor = self.o.coor()
      self.assert_list_almost_equal(result_coor[0][10],expected_coor_sample,3)


   def test_1AA_NoEND(self):
      '''
	   test a 1AA pdb file with 1frame and without END statement
	   '''
      #
      self.o.read_pdb(moduleDataPath+'1AA-NoEND.pdb')
      result_coor = self.o.coor()
      result_sum_coor = sum(sum(sum((result_coor))))
      #
      expected_coor = numpy.array([[-21.525, -67.562,  86.759], [-22.003, -68.460, 86.892],[-21.905, -66.929,  87.525],[-20.492, -67.726, 86.876],[-21.725, -66.910, 85.457],[-21.476, -67.600, 84.661],[-21.157, -65.997, 85.450],[-23.103, -66.411, 85.215],[-23.249, -65.504, 84.385]],floattype)
      expected_sum_coor = sum(sum(expected_coor))
      #
      self.assertAlmostEqual(result_sum_coor,expected_sum_coor,3)

   def test_cleaned_up_package_rna(self):
      '''
	   test a pdb file of rna with 1250 frames of size 1.0g
      '''
      #
      self.o.read_pdb(moduleDataPath+"new_package_rna.pdb")
      result_coor = self.o.coor()
      print '\nlength of result_coor \n',len(result_coor[0])
      print '\nresult_coor \n',result_coor
      #
      self.assertEqual(len(result_coor[0]),3719)
      expected_coor_sample = numpy.array([-12.872, 13.360, -153.873],floattype) #atom 10627 of frame 3
      self.assert_list_almost_equal(result_coor[0][299],expected_coor_sample,3)

   def test_problem_pdb(self):
      '''
	   test a pdb file with non-charmm atom names
      '''
      #
      print 'ZHL'
      self.o.read_pdb(moduleDataPath+"nef_nohis.pdb")
      print self.o.name()

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

