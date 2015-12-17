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

class Test_intg_sasmol_SasSys_remove_object(MockerTestCase):

   def setUp(self):
      pass

   def test_null_remove_nothing_from_atomic(self):
      '''
      null test remove nothing from the atomic type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.remove_object()
      self.assertEqual(o._objectarray[0][1].natoms(),0)

   def test_remove_atomic_from_atomic(self):
      '''
      test remove atomic from the atomic type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.remove_object(id=0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)

   def test_remove_sol_from_sol(self):
      '''
      test remove solid from the solid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('solid')
      o.remove_object(id=0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1]._name,'Sol_None')

   def test_remove_hybrid_from_hybrid(self):
      '''
      test remove hybrid from the hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('hybrid')
      o.remove_object(id=0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1]._name,'Sol_None')

   def test_remove_atomic_from_atomic_sol_hybrid(self):
      '''
      test remove atomic from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_sol_from_atomic_sol_hybrid(self):
      '''
      test remove solid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=1)
      self.assertEqual(o._objectarray[0][1].natoms(),0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')


   def test_remove_hybrid_from_atomic_sol_hybrid(self):
      '''
      test remove hybrid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=2)
      self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_atomic_from_atomic_sol_hybrid(self):
      '''
      test remove atomic from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_atomic_sol_from_atomic_sol_hybrid(self):
      '''
      test remove atomic_solid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=0)
      o.remove_object(id=1)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_atomic_hybrid_from_atomic_sol_hybrid(self):
      '''
      test remove atomic_solid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=0)
      o.remove_object(id=2)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)
      self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_solid_hybrid_from_atomic_sol_hybrid(self):
      '''
      test remove solid_hybrid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=1)
      o.remove_object(id=2)
      self.assertEqual(o._objectarray[0][1].natoms(),0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_remove_atomic_solid_hybrid_from_atomic_sol_hybrid(self):
      '''
      test remove atomic_solid_hybrid from the atomic+solid_hybrid type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.add_object('solid')
      o.add_object('hybrid')
      o.remove_object(id=0)
      o.remove_object(id=1)
      o.remove_object(id=2)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[0][1].natoms(),0)
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[1][1]._name,'Sol_None')
      with self.assertRaises(Exception): self.assertEqual(o._objectarray[2][1]._name,'Hybrid_None')

   def test_negative_remove_sol_from_atomic(self):
      '''
      negative test remove solid from the atomic type
      '''
      #
      o=sasmol.SasSys(0)
      o.add_object('atomic')
      o.remove_object(id=1)
      self.assertEqual(o._objectarray[0][1].natoms(),0)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

