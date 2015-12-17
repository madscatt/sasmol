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

from unittest import main 
from mocker import Mocker, MockerTestCase
import sasmol.sasmol as sasmol

import os
DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasmol')+os.path.sep

class Test_intg_sasmol_SasAtm_name(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasAtm(3,'1CRN-3frames.pdb')

   def test_1CRN_3frames(self):
      '''
	   test a regular pdb file with 3 frame
	   '''
      #
      expected = ['N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'N', 'CA', 'C', 'O', 'CB', 'OG', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'N', 'CA', 'C', 'O', 'CB', 'OG', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'NE', 'CZ', 'NH1', 'NH2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'N', 'CA', 'C', 'O', 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'OE1', 'OE2', 'N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'CA', 'C', 'O', 'CB', 'CG1', 'CG2', 'CD1', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'N', 'CA', 'C', 'O', 'N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB', 'OG1', 'CG2', 'N', 'CA', 'C', 'O', 'CB', 'SG', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD', 'N', 'CA', 'C', 'O', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'OD2', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ', 'OH', 'N', 'CA', 'C', 'O', 'CB', 'N', 'CA', 'C', 'O', 'CB', 'CG', 'OD1', 'ND2', 'OXT']
      #
      self.o.read_pdb(DataPath+'1CRN-3frames.pdb')
      #
      result = self.o.name()
      print result
      #
      self.assertEqual(expected, result)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

