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

'''
Contract for unit test of sasio_Files_element_filter:

Test for null
Test for conflict atoms in the case of heavy atoms
Test for conflict atoms in the case of light atoms
Test for H/C/N/O/S/P atoms
Test for miscellaneous atoms (CAL, POT, ...)
Test for noncharmm/wrong atoms (ABC, ...)
'''


from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import os

# This data for atomic properties are stored under sasproperties folder

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasproperties')+os.path.sep

class Test_unit_sasio_Files_get_elements(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)


   def test_null(self):
      '''
	   test for null atoms
      '''
      #
      name = '  '
      resname = ''
      (error, element_name) = self.o.get_element(name,resname)
      self.assertTrue(len(error)>0)
      self.assertEqual(element_name,'')

   def test_conflict_atoms_heavy_element(self):
      '''
	   test for conflict atoms when they are heavy elments
      '''
      #
      datafile = DataPath+'ConflictAtoms.txt'
      resname = 'RES'
      ltmp = []
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         resname = name
         realname = name
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)


   def test_conflict_atoms_light_element(self):
      '''
	   test for conflict atoms when they are light elments
      '''
      #
      datafile = DataPath+'ConflictAtoms.txt'
      resname = 'RES'
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         realname = atom.split()[1]
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)


   def test_H(self):
      '''
      test for H atoms
      '''
      #
      datafile = DataPath+'Hatoms.txt'
      resname = 'RES'
      realname = 'H'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_C(self):
      '''
      test for C atoms
      '''
      #
      datafile = DataPath+'Catoms.txt'
      resname = 'RES'
      realname = 'C'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_N(self):
      '''
      test for N atoms
      '''
      #
      datafile = DataPath+'Natoms.txt'
      resname = 'RES'
      realname = 'N'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_O(self):
      '''
      test for O atoms
      '''
      #
      datafile = DataPath+'Oatoms.txt'
      resname = 'RES'
      realname = 'O'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_S(self):
      '''
      test for S atoms
      '''
      #
      datafile = DataPath+'Satoms.txt'
      resname = 'RES'
      realname = 'S'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_P(self):
      '''
      test for P atoms
      '''
      #
      datafile = DataPath+'Patoms.txt'
      resname = 'RES'
      realname = 'P'
      for atom in open(datafile).readlines():
         name = atom.strip()
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_Other(self):
      '''
      test for periodic table heavy atoms
      '''
      #
      datafile = DataPath+'Otheratoms_2.txt'
      resname = 'RES'
      for atom in open(datafile).readlines():
         name = atom.strip()
         realname = name
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name,realname
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)


   def test_miscellaneous(self):
      '''
      test for miscellaneous atoms
      '''
      #
      datafile = DataPath+'MisAtoms.txt'
      resname = 'RES'
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         realname = atom.split()[1]
         (error, element_name) = self.o.get_element(name,resname)
         print error, element_name
         self.assertEqual(len(error),0)
         self.assertEqual(element_name,realname)

   def test_non_charmm1(self):
      '''
      test for non-charmm atoms-OXT
      '''
      #
      name = 'OXT'
      resname = 'RES'
      (error, element_name) = self.o.get_element(name,resname)
      print error, element_name
      self.assertTrue(len(error)==0)
      self.assertEqual(element_name,'O')

   def test_non_charmm2(self):
      '''
      test for non-charmm atoms-1H11
      '''
      #
      name = '1H11'
      resname = 'RES'
      (error, element_name) = self.o.get_element(name,resname)
      print error, element_name
      self.assertTrue(len(error)==0)
      self.assertEqual(element_name,'H')


   def test_wrong1(self):
      '''
      test for non-charmm or wrong atoms
      '''
      #
      name = '$ABC'
      resname = 'RES'
      (error, element_name) = self.o.get_element(name,resname)
      print error, element_name
      self.assertTrue(len(error)>0)
      self.assertEqual(element_name,'')

   def test_wrong1(self):
      '''
      test for non-charmm or wrong atoms
      '''
      #
      name = '11H'
      resname = 'RES'
      (error, element_name) = self.o.get_element(name,resname)
      print error, element_name
      self.assertTrue(len(error)>0)
      self.assertEqual(element_name,'')


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

