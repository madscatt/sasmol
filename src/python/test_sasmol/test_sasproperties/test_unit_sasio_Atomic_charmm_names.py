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


contract:

make sure the list is unique
make sure the right hydrogen list was generated
make sure the right carbon list was generated
make sure the right nitrogen list was generated
make sure the right oxygen list was generated
make sure the right sulfur list was generated
make sure the right phosphorous list was generated
make sure the right other list (from the periodic table) was generated
'''

from unittest import main 
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasproperties')+os.path.sep


class Test_unit_sasproperties_Atomic_charmm_names(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)
      (self.hl, self.cl, self.nl, self.ol, self.sl, self.pl, self.otherl)=self.o.charmm_names()

   def unique(self,seq):
      '''
      remove the duplicate elements
      '''
      seen = set()
      seen_add = seen.add
      return [ x for x in seq if x not in seen and not seen_add(x)]

   def compare_namelists(self, l1, l2):
      '''
      check whether two list contains the same atoms
      the order doesnt matter
      '''
      if len(l1)!=len(l2):
         return -1
      else:
         for i1 in l1:
            if i1 not in l2:
               return -1
      return 0


   def test_uniqe(self):
      '''
	   make sure the list is unique
      '''
      #
      lall = self.hl + self.cl + self.nl + self.ol + self.sl + self.pl + self.otherl
      lallunique = self.unique(lall)
      self.assertEqual(lall,lallunique)



   def test_H(self):
      '''
      make sure the right hydrogen list was generated
      '''
      #
      datafile = DataPath+'Hatoms.txt'
      hltmp = []
      for atom in open(datafile).readlines():
         hltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(hltmp, self.hl), 0)


   def test_C(self):
      '''
      make sure the right carbon list was generated
      '''
      #
      datafile = DataPath+'Catoms.txt'
      cltmp = []
      for atom in open(datafile).readlines():
         cltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(cltmp, self.cl), 0)


   def test_N(self):
      '''
      make sure the right nitrogen list was generated
      '''
      #
      datafile = DataPath+'Natoms.txt'
      nltmp = []
      for atom in open(datafile).readlines():
         nltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(nltmp, self.nl), 0)


   def test_O(self):
      '''
      make sure the right oxygen list was generated
      '''
      #
      datafile = DataPath+'Oatoms.txt'
      oltmp = []
      for atom in open(datafile).readlines():
         oltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(oltmp, self.ol), 0)

   def test_S(self):
      '''
      make sure the right sulfur list was generated
      '''
      #
      datafile = DataPath+'Satoms.txt'
      sltmp = []
      for atom in open(datafile).readlines():
         sltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(sltmp, self.sl), 0)


   def test_P(self):
      '''
      make sure the right phosphorous list was generated
      '''
      #
      datafile = DataPath+'Patoms.txt'
      pltmp = []
      for atom in open(datafile).readlines():
         pltmp.append(atom.strip())
      self.assertEqual(self.compare_namelists(pltmp, self.pl), 0)


   def test_Other(self):
      '''
      make sure the right other list (from the periodic table) was generated
      '''
      #
      datafile = DataPath+'Otheratoms.txt'
      otherltmp = []
      for atom in open(datafile).readlines():
         otherltmp.append(atom.strip())
      #print 'otherltmp = ',otherltmp
      #print 'self.otherl = ',self.otherl
      #print list(set(self.otherl) - set(otherltmp))
      self.assertEqual(self.compare_namelists(otherltmp, self.otherl), 0)


   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

