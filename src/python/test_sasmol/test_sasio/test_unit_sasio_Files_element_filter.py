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
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol

import os, sys, string

# This data for atomic properties are stored under sasproperties folder

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasproperties')+os.path.sep

class Test_unit_sasio_Files_element_filter(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)


   def test_null(self):
      '''
	   test for null atom list
      '''
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      self.o.element_filter()
      result_elements = self.o.element()
      self.assertTrue(expected_elements==result_elements)


   def test_conflict_atoms_heavy_element(self):
      '''
	   test for conflict atoms when they are heavy elments
      '''
      #
      datafile = DataPath+'ConflictAtoms.txt'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         resname = name
         element = name
         names.append(name)
         resnames.append(resname)
         elements.append('  ')
         expected_elements.append(element)
      print 'names\n',names
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      self.o.element_filter()
      result_elements = self.o.element()
      self.assertTrue(expected_elements==result_elements)


   def test_conflict_atoms_light_element(self):
      '''
	   test for conflict atoms when they are light elments
      '''
      #
      datafile = DataPath+'ConflictAtoms.txt'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         element = atom.split()[1]
         resname = 'RES'
         names.append(name)
         resnames.append(resname)
         elements.append('  ')
         expected_elements.append(element)
      print 'names\n',names
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      self.o.element_filter()
      result_elements = self.o.element()
      self.assertTrue(expected_elements==result_elements)


   def test_H(self):
      '''
      test for H atoms
      '''
      #
      datafile = DataPath+'Hatoms.txt'
      testelement = 'H'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()


   def test_C(self):
      '''
      test for C atoms
      '''
      #
      datafile = DataPath+'Catoms.txt'
      testelement = 'C'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()

   
   def test_N(self):
      '''
      test for N atoms
      '''
      #
      datafile = DataPath+'Natoms.txt'
      testelement = 'N'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()



   def test_O(self):
      '''
      test for O atoms
      '''
      #
      datafile = DataPath+'Oatoms.txt'
      testelement = 'O'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()


   def test_S(self):
      '''
      test for S atoms
      '''
      #
      datafile = DataPath+'Satoms.txt'
      testelement = 'S'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()


   def test_P(self):
      '''
      test for P atoms
      '''
      #
      datafile = DataPath+'Patoms.txt'
      testelement = 'P'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      if len(names)>0:
         self.o.setName(names)
         self.o.setResname(resnames)
         self.o.setElement(elements)
         self.o.element_filter()
         result_elements = self.o.element()
         self.assertTrue(expected_elements==result_elements)
      #
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         if not name.isupper():
            element = testelement
            resname = 'RES'
            names.append(name)
            resnames.append(resname)
            elements.append('  ')
            expected_elements.append(element) 
      print 'names\n',names
      self.o.setName(names)
      if len(names)>0:
         self.o.setResname(resnames)
         self.o.setElement(elements)
         with self.assertRaises(SystemExit):
            self.o.element_filter()


   """
   def test_Other(self):
      '''
      test for periodic table heavy atoms
      '''
      #
      datafile = DataPath+'Otheratoms.txt'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.strip()
         element = name
         resname = 'RES'
         names.append(name)
         resnames.append(resname)
         elements.append('  ')
         expected_elements.append(element)
      print 'names\n',names
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      self.o.element_filter()
      result_elements = self.o.element()
      print '\nexpected_elements:\n',expected_elements
      print '\nresult_elements:\n',result_elements
      self.assertTrue(expected_elements==result_elements)
   """



   def test_miscellaneous(self):
      '''
      test for miscellaneous atoms
      '''
      #
      datafile = DataPath+'MisAtoms.txt'
      names = []
      resnames = []
      elements = []
      expected_elements = []
      for atom in open(datafile).readlines():
         name = atom.split()[0]
         element = atom.split()[1]
         resname = 'RES'
         names.append(name)
         resnames.append(resname)
         elements.append('  ')
         expected_elements.append(element)
      print 'names\n',names
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      self.o.element_filter()
      result_elements = self.o.element()
      print '\nexpected_elements:\n',expected_elements
      print '\nresult_elements:\n',result_elements
      self.assertTrue(expected_elements==result_elements)

   
   def test_wrong(self):
      '''
      test for non-charmm or wrong atoms
      '''
      #
      names = ['H','ABC','C']
      resnames = ['RES','RES','RES']
      elements = [' ',' ',' ']
      self.o.setName(names)
      self.o.setResname(resnames)
      self.o.setElement(elements)
      #
      self.o.element_filter()


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

