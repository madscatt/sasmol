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


class Test_intg_sasio_Files_initialize_children(MockerTestCase):

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

   def test_null(self):
      '''
      null test by not reading pdb
      '''
      #
      with self.assertRaises(Exception):
         self.o.initialize_children()


   def test_1ATM(self):
      '''
      test everything for a pdb file with 1 atom and 1 frame
      '''
      #
      self.o.read_pdb(DataPath+'1ATM.pdb')
      self.o.initialize_children()
      self.assert_list_almost_equal(self.o.names_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.resnames_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.resids_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.chains_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.occupancies_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.betas_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.elements_mask(), numpy.array([[1]],numpy.int),3)
      self.assert_list_almost_equal(self.o.segnames_mask(), numpy.array([[1]],numpy.int),3)

   def test_2AAD_names(self):
      '''
      test names for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_names_mask = self.o.names_mask()
      print '\nresult_names_mask \n',result_names_mask.tolist()
      #
      expected_names_mask = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0]],numpy.int)

      print '\nexpected_names_mask \n',expected_names_mask
      #
      self.assert_list_almost_equal(expected_names_mask, result_names_mask,3)
 
   def test_2AAD_resnames(self):
      '''
      test resnames for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.resnames_mask()
      print '\nresult_names_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]],numpy.int)

      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_resids(self):
      '''
      test resids for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.resids_mask()
      print '\nresult_names_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]],numpy.int)

      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_chains(self):
      '''
      test chains for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.chains_mask()
      print '\nresult_names_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],numpy.int)
      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_occupancies(self):
      '''
      test occupancies for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.occupancies_mask()
      print '\nresult_names_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],numpy.int)
      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_betas(self):
      '''
      test betas for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.betas_mask()
      print '\nresult_names_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0], [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1]],numpy.int)
      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_elements(self):
      '''
      test elements for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.elements_mask()
      print '\nresult_elements_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0], [0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 1], [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0]],numpy.int)
      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   def test_2AAD_senames(self):
      '''
      test segnames for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_mask = self.o.segnames_mask()
      print '\nresult_segnames_mask \n',result_mask.tolist()
      #
      expected_mask = numpy.array([[1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]],numpy.int)
      print '\nexpected_mask \n',expected_mask
      #
      self.assert_list_almost_equal(expected_mask, result_mask,3)

   
   def test_rna(self):
      '''
      test a pdb file of rna
      '''
      #
      self.o.read_pdb(DataPath+"rna.pdb")
      self.o.initialize_children()
      self.assertEqual(sum(sum(self.o.names_mask())), 10632)
      self.assertEqual(sum(sum(self.o.resnames_mask())), 10632)
      self.assertEqual(sum(sum(self.o.resids_mask())), 10632)
      self.assertEqual(sum(sum(self.o.occupancies_mask())), 10632)
      self.assertEqual(sum(sum(self.o.chains_mask())), 10632)
      self.assertEqual(sum(sum(self.o.betas_mask())), 10632)
      self.assertEqual(sum(sum(self.o.elements_mask())), 10632)
      self.assertEqual(sum(sum(self.o.segnames_mask())), 10632)

   def test_1CRN(self):
      '''
      test a small protein (1CRN)
      '''
      #
      self.o.read_pdb(DataPath+"1CRN.pdb")
      self.o.initialize_children()
      self.assertEqual(sum(sum(self.o.names_mask())), 327)
      self.assertEqual(sum(sum(self.o.resnames_mask())), 327)
      self.assertEqual(sum(sum(self.o.resids_mask())), 327)
      self.assertEqual(sum(sum(self.o.occupancies_mask())), 327)
      self.assertEqual(sum(sum(self.o.chains_mask())), 327)
      self.assertEqual(sum(sum(self.o.betas_mask())), 327)
      self.assertEqual(sum(sum(self.o.elements_mask())), 327)
      self.assertEqual(sum(sum(self.o.segnames_mask())), 327)

   #@skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")   
   def test_1KP8(self):
      '''
      test a pdb file of 1KP8
      '''
      #
      self.o.read_pdb(DataPath+"1KP8.pdb")
     # with self.assertRaises(Exception):
     #     self.o.initialize_children()
      self.o.initialize_children()
      
      self.assertEqual(sum(sum(self.o.names_mask())), 57085)
      self.assertEqual(sum(sum(self.o.resnames_mask())), 57085)
      self.assertEqual(sum(sum(self.o.resids_mask())), 57085)
      self.assertEqual(sum(sum(self.o.occupancies_mask())), 57085)
      self.assertEqual(sum(sum(self.o.chains_mask())), 57085)
      self.assertEqual(sum(sum(self.o.betas_mask())), 57085)
      self.assertEqual(sum(sum(self.o.elements_mask())), 57085)
      self.assertEqual(sum(sum(self.o.segnames_mask())), 57085)
     

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

