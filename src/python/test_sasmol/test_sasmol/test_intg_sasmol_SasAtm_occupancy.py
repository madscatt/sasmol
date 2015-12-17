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
import sasmol.sasmol as sasmol
from unittest import main 
from mocker import Mocker, MockerTestCase

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasmol')+os.path.sep

class Test_intg_sasmol_SasAtm_occupancy(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasAtm(3,'1CRN-3frames.pdb')

   def test_1CRN_3frames(self):
      '''
	   test a regular pdb file with 3 frame
	   '''
      #
      expected = ['1.00']*327
      #
      self.o.read_pdb(DataPath+'1CRN-3frames.pdb')
      #
      result = self.o.occupancy()
      print result
      #
      self.assertEqual(expected, result)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

