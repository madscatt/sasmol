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

class Test_unit_sasio_Files_print_error(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)

   def test_no_error(self):
      '''
      test for no error
      '''
      error = []
      self.o.check_error(error)

   def test_error(self):
      '''
      test for error
      '''
      error = ['wrong']
      with self.assertRaises(SystemExit):
         self.o.check_error(error)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

