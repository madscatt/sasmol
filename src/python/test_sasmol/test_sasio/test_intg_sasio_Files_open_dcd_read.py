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

from unittest import main, skipIf
from mocker import Mocker, MockerTestCase

import sasmol.sasmol as sasmol

import os

DataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep

class Test_intg_sasio_Files_open_dcd_read(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)



   def test_1ATM(self):
      '''
	   test a dcd with 2 frames based on a 1-atom pdb
	   '''
      #
      dcdFileName = DataPath+'1ATM.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],1)
      self.assertEqual(fp[2],2)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],5)

   def test_2AAD(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb
	   '''
      #
      dcdFileName = DataPath+'2AAD.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],15)
      self.assertEqual(fp[2],3)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],5)


   def test_rna_1to10frames(self):
      '''
	   test a dcd of 10 frames based on a rna molecule
	   '''
      #
      dcdFileName = DataPath+'rna-1to10.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],10632)
      self.assertEqual(fp[2],10)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],5)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
   def test_rna_1point0g(self):
      '''
	   test a dcd of size 1.0gb based on a rna molecule
	   '''
      #
      dcdFileName = '/tmp/rna-1.0g.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],10632)
      self.assertEqual(fp[2],7813)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],0)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
   def test_rna_2point0g(self):
      '''
	   test a dcd of size 2.0gb based on a rna molecule
	   '''
      #
      dcdFileName = '/tmp/rna-2.0g.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],10632)
      self.assertEqual(fp[2],15625)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],0)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
   def test_rna_3point2g(self):
      '''
	   test a dcd of size 3.2gb based on a rna molecule
	   '''
      #
      dcdFileName = '/tmp/rna-3.2g.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],10632)
      self.assertEqual(fp[2],25000)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],0)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")
   def test_rna_6point4g(self):
      '''
	   test a dcd of size 6.4gb based on a rna molecule
	   '''
      #
      dcdFileName = '/tmp/rna-6.4g.dcd'
      fp = self.o.open_dcd_read(dcdFileName)
      self.assertEqual(str(type(fp[0])),"<type 'SwigPyObject'>")
      self.assertEqual(fp[1],10632)
      self.assertEqual(fp[2],50000)
      self.assertEqual(fp[3],0)
      self.assertEqual(fp[4],0)

   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

