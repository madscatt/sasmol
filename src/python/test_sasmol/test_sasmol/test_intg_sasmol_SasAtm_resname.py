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

class Test_intg_sasmol_SasAtm_resname(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasAtm(3,'1CRN-3frames.pdb')

   def test_1CRN_3frames(self):
      '''
	   test a regular pdb file with 3 frame
	   '''
      #
      expected = ['THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'SER', 'SER', 'SER', 'SER', 'SER', 'SER', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'SER', 'SER', 'SER', 'SER', 'SER', 'SER', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'PHE', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'VAL', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'ARG', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'LEU', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'GLY', 'GLY', 'GLY', 'GLY', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'GLU', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'GLY', 'GLY', 'GLY', 'GLY', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'ILE', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'GLY', 'GLY', 'GLY', 'GLY', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'THR', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'CYS', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'PRO', 'GLY', 'GLY', 'GLY', 'GLY', 'ASP', 'ASP', 'ASP', 'ASP', 'ASP', 'ASP', 'ASP', 'ASP', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'TYR', 'ALA', 'ALA', 'ALA', 'ALA', 'ALA', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN', 'ASN']
      #
      self.o.read_pdb(DataPath+'1CRN-3frames.pdb')
      #
      result = self.o.resname()
      print result
      #
      self.assertEqual(expected, result)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

