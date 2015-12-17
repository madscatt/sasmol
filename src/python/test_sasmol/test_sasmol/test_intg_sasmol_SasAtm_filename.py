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
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol
import sasmol.sasop as sasop
import sasmol.sascalc as sascalc

import numpy, os, copy

import warnings; warnings.filterwarnings('ignore')

floattype=os.environ['SASSIE_FLOATTYPE']

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasmol')+os.path.sep

class Test_intg_sasmol_SasAtm_Filename(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasAtm(3,'1CRN-3frames.pdb')

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

   def test_1CRN_3frames(self):
      '''
	   test a regular pdb file with 3 frame
	   '''
      #
      fileName = '1CRN-3frames.pdb'
      expected_fileName = fileName
      #
      self.o.read_pdb(DataPath+'1CRN-3frames.pdb')
      self.o.setFilename(fileName)
      #
      result_fileName = self.o.filename()
      #
      self.assertEqual(expected_fileName, result_fileName)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

