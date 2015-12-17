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

class Test_intg_sasmol_SasSys_add_object(MockerTestCase):

   def setUp(self):
      pass


   def test_id_systype_atomic(self):
      '''
      test add object for an atomic type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      self.assertEqual(o._objectarray[0][0],0)
      self.assertEqual(o._objectarray[0][1].id(),0)
      self.assertEqual(o._objectarray[0][1].totalmass(),0.0)
      self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[0][1].mass(),None)
      self.assertEqual(o._objectarray[0][1].coor(),None)
      self.assertEqual(o._objectarray[0][1].com(),None)

   def test_id_systype_sol(self):
      '''
      test add object for a solid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('solid')
      self.assertEqual(o._objectarray[0][0],0)
      self.assertEqual(o._objectarray[0][1]._name,'Sol_None')

   def test_id_systype_hybrid(self):
      '''
      test add object for a hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('hybrid')
      self.assertEqual(o._objectarray[0][0],0)
      self.assertEqual(o._objectarray[0][1]._name,'Hybrid_None')

   def test_id_systype_atomic_sol_hybrid(self):
      '''
      test add object for atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      self.assertEqual(o._objectarray[0][0],0)
      self.assertEqual(o._objectarray[0][1].id(),0)
      self.assertEqual(o._objectarray[0][1].totalmass(),0.0)
      self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[0][1].mass(),None)
      self.assertEqual(o._objectarray[0][1].coor(),None)
      self.assertEqual(o._objectarray[0][1].com(),None)
      self.assertEqual(o._objectarray[1][0],1)
      self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      self.assertEqual(o._objectarray[2][0],2)
      self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_wrong(self):
      '''
      test add object for a wrong type
      '''
      #
      import sys,os
      o=sasmol.SasSys(0)
      stdoutFileName = __file__+'.stdiout'
      sys.stdout = open(stdoutFileName,'w')
      o.add_object('whatisthis')
      sys.stdout = sys.__stdout__
      info=open(stdoutFileName).readlines()
      expected_info = ['>>> error: need to specify addtype == atom, molecule, or assembly\n','>>> no objects created\n']
      self.assertEqual(info,expected_info)
      os.remove(stdoutFileName)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

