'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 
	 Core-Testing: Copyright (C) 2011 Hailiang Zhang, Ph.D.

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


from unittest import main 
from mocker import Mocker, MockerTestCase
from sassie.sasmol import sasmol





class Test_intg_sasmol_SasSys_init(MockerTestCase):

   def setUp(self):
      pass


   def test_id(self):
      '''
      test initializer with id
      '''
      #
      id=3
      o=sasmol.SasSys(id)
      self.assertEqual(o.id(),id)
      self.assertEqual(o.systype(),'atomic')


   def test_id_systype_atomic(self):
      '''
      test initializer with id and systype of atomic
      '''
      #
      id=3
      systype='atomic'
      o=sasmol.SasSys(id,systype)
      self.assertEqual(o.id(),id)
      self.assertEqual(o.systype(),systype)

   def test_id_systype_solid(self):
      '''
      test initializer with id and systype of solid
      '''
      #
      id=3
      systype='solid'
      o=sasmol.SasSys(id,systype)
      self.assertEqual(o.id(),id)
      self.assertEqual(o.systype(),systype)

   """
   def test_id_systype_hybrid(self):
      '''
      test initializer with id and systype of solid
      '''
      #
      id=3
      systype='hybrid'
      o=sasmol.SasSys(id,systype)
      self.assertEqual(o.id(),id)
      self.assertEqual(o.systype(),systype)
   """

   def test_id_wrong(self):
      '''
      test initializer with the wrong systype input
      '''
      #
      id=3
      systype='whatisthis'
      o=sasmol.SasSys(id,systype)
      self.assertEqual(o.id(),id)
      with self.assertRaises(Exception):
         self.assertEqual(o.systype(),systype)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 
