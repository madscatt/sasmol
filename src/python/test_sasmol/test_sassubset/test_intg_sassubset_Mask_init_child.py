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
      null test by providing the empty parameter to init_child
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      with self.assertRaises(Exception):
         self.o.init_child('')

   def test_negative(self):
      '''
      negative test by providing the wrong parameter to init_child
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      with self.assertRaises(Exception):
         self.o.init_child('ab')

   def test_1ATM(self):
      '''
      test everything for a pdb file with 1 atom and 1 frame
      '''
      #
      self.o.read_pdb(DataPath+'1ATM.pdb')
      self.o.initialize_children()
      self.assertEqual([item.name() for item in self.o.init_child('names')],[['N']])
      self.assertEqual([item.resname() for item in self.o.init_child('resnames')],[['ILE']])
      self.assertEqual([item.resid() for item in self.o.init_child('resids')],[[515]])
      self.assertEqual([item.chain() for item in self.o.init_child('chains')],[['N']])
      self.assertEqual([item.occupancy() for item in self.o.init_child('occupancies')],[['1.00']])
      self.assertEqual([item.beta() for item in self.o.init_child('betas')],[['36.37']])
      #self.assertEqual([item.segname() for item in self.o.init_child('segnames')],[['DUM0']])
      self.assertEqual([item.segname() for item in self.o.init_child('segnames')],[['N']])
      self.assertEqual([item.element() for item in self.o.init_child('elements')],[['N']])

   def test_2AAD_names(self):
      '''
      test names for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('names')
      result = []
      for item in result_child:
          result.append(item.name())
      print 'result \n',result
      #
      expected = [['N', 'N'] , ['CA', 'CA'] , ['C', 'C'] , ['O', 'O'] , ['CB', 'CB'] , ['CG1'] , ['CG2', 'CG2'] , ['CD1'] , ['OG1']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)
 
   def test_2AAD_resnames(self):
      '''
      test resnames for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('resnames')
      result = []
      for item in result_child:
          result.append(item.resname())
      print 'result \n',result
      #
      expected = [['ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE'], ['THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_resids(self):
      '''
      test resids for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('resids')
      result = []
      for item in result_child:
          result.append(item.resid().tolist())
      print 'result \n',result
      #
      expected = [[515, 515, 515, 515, 515, 515, 515, 515],[516, 516, 516, 516, 516, 516, 516]]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_chains(self):
      '''
      test chains for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('chains')
      result = []
      for item in result_child:
          result.append(item.chain())
      print 'result \n',result
      #
      expected = [['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_occupancies(self):
      '''
      test occupancies for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('occupancies')
      result = []
      for item in result_child:
          result.append(item.occupancy())
      print 'result \n',result
      #
      expected = [['1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00', '1.00']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_betas(self):
      '''
      test betas for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('betas')
      result = []
      for item in result_child:
          result.append(item.beta())
      print 'result \n',result
      #
      expected = [['36.37'], ['36.23'], ['36.32'], ['36.04'], ['36.69'], ['38.12'], ['34.42'], ['39.85'], ['35.01'], ['35.51'], ['38.09'], ['36.94'], ['36.74'], ['37.19'], ['34.44']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_elements(self):
      '''
      test elements for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('elements')
      result = []
      for item in result_child:
          result.append(item.element())
      print 'result \n',result
      #
      expected = [['N', 'N'], ['C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C'], ['O', 'O', 'O']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_senames(self):
      '''
      test segnames for a pdb file with 2 amino acids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('segnames')
      result = []
      for item in result_child:
          result.append(item.segname())
      print 'result \n',result
      #
      #expected = [['DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0', 'DUM0']]
      expected = [['N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   def test_2AAD_mix(self):
      '''
      test for elements for a pdb file with 2 amino acids which has child grouped by resids
      '''
      #
      self.o.read_pdb(DataPath+'2AAD.pdb')
      self.o.initialize_children()
      result_child = self.o.init_child('resids')
      result = []
      for item in result_child:
          result.append(item.element())
      print 'result \n',result
      #
      expected = [['N', 'C', 'C', 'O', 'C', 'C', 'C', 'C'], ['N', 'C', 'C', 'O', 'C', 'O', 'C']]
      print 'expected \n',expected
      #
      self.assertEqual(expected, result)

   
   def test_rna(self):
      '''
      test a pdb file of rna
      '''
      #
      self.o.read_pdb(DataPath+"rna.pdb")
      self.o.initialize_children()
      result_child = self.o.init_child('resids')
      self.assertEqual(len(result_child),25)
      sample=2
      self.assertEqual(result_child[sample].resid().tolist(),[3]*520)
      self.assertEqual(sum([sum(item.resid()) for item in result_child]),118880)

   def test_1CRN(self):
      '''
      test a small protein (1CRN)
      '''
      #
      self.o.read_pdb(DataPath+"1CRN.pdb")
      self.o.initialize_children()
      result_child = self.o.init_child('resids')
      self.assertEqual(len(result_child),46)
      sample=2
      self.assertEqual(result_child[sample].resid().tolist(),[3]*6)
      self.assertEqual(sum([sum(item.resid()) for item in result_child]),7637)

   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")   
   def test_1KP8(self):
      '''
      test a pdb file of 1KP8
      '''
      #
      self.o.read_pdb(DataPath+"1KP8.pdb")
      with self.assertRaises(Exception):
          self.o.initialize_children()
      '''
      self.assertEqual(sum(sum(self.o.names_mask())), 57085)
      self.assertEqual(sum(sum(self.o.resnames_mask())), 57085)
      self.assertEqual(sum(sum(self.o.resids_mask())), 57085)
      self.assertEqual(sum(sum(self.o.occupancies_mask())), 57085)
      self.assertEqual(sum(sum(self.o.chains_mask())), 57085)
      self.assertEqual(sum(sum(self.o.betas_mask())), 57085)
      self.assertEqual(sum(sum(self.o.elements_mask())), 57085)
      self.assertEqual(sum(sum(self.o.segnames_mask())), 57085)
      '''

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

