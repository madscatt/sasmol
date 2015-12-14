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
import sasmol.dcdio as dcdio

import os, sys, string, shutil

pdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep
dcdDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','dcd_common')+os.path.sep
moduleDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','sasmol','sasio')+os.path.sep

class Test_intg_sasio_Files_open_dcd_write(MockerTestCase):

   def setUp(self):
      self.o=sasmol.SasMol(0)


   def test_file_doesnt_exist(self):
      '''
	   test a dcd which doent exist
	   '''
      filename = 'file-notexist.dcd'
      pdbFile = pdbDataPath+'1ATM.pdb'
      dcdFile = moduleDataPath+'test-results/'+filename
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(dcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(dcdFile)


   def test_1ATM(self):
      '''
	   test a dcd with 2 frames based on an 1-atom pdb
	   '''
      #
      filecode = '1ATM'
      pdbFile = pdbDataPath+filecode+'.pdb'
      dcdFile = dcdDataPath+filecode+'.dcd'
      tmpDcdFile = moduleDataPath+'test-results/'+filecode+'.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   def test_2AAD(self):
      '''
	   test a dcd with 3 frames based on a 2-aa pdb

	   '''
      #
      filecode = '2AAD'
      pdbFile = pdbDataPath+filecode+'.pdb'
      dcdFile = dcdDataPath+filecode+'.dcd'
      tmpDcdFile = moduleDataPath+'test-results/'+filecode+'.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   def test_rna_1to10(self):
      '''
	   test a dcd with 10 frames based on a 2-aa pdb

	   '''
      #
      filecode = 'rna-1to10'
      pdbFile = pdbDataPath+filecode+'.pdb'
      dcdFile = dcdDataPath+filecode+'.dcd'
      tmpDcdFile = moduleDataPath+'test-results/'+filecode+'.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_0point8gb(self):
      '''
	   test a dcd of size 0.8gb based on a rna molecule
	   '''
      #
      filename = "rna-0.8g.dcd"
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-0.8g.dcd'
      tmpDcdFile = moduleDataPath+'test-results/rna-0.8g.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_1point0gb(self):
      '''
	   test a dcd of size 1.0gb based on a rna molecule
	   '''
      #
      filename = "rna-1.0g.dcd"
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-1.0g.dcd'
      tmpDcdFile = moduleDataPath+'test-results/rna-1.0g.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)



   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_2point0gb(self):
      '''
	   test a dcd of size 2.0gb based on a rna molecule
	   '''
      #
      filename = "rna-2.0g.dcd"
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-2.0g.dcd'
      tmpDcdFile = moduleDataPath+'test-results/rna-2.0g.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   
   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_3point2gb(self):
      '''
	   test a dcd of size 3.2gb based on a rna molecule
	   '''
      #
      filename = "rna-3.2g.dcd"
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-3.2g.dcd'
      tmpDcdFile = moduleDataPath+'test-results/rna-3.2g.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)


   
   @skipIf(os.environ['SASSIE_HUGETEST']=='n',"I am not testing huge files")   
   def test_rna_6point4gb(self):
      '''
	   test a dcd of size 6.4gb based on a rna molecule
	   '''
      #
      filename = "rna-6.4g.dcd"
      pdbFile = pdbDataPath+'rna.pdb'
      dcdFile = '/tmp/rna-6.4g.dcd'
      tmpDcdFile = moduleDataPath+'test-results/rna-6.4g.dcd'
      shutil.copy(dcdFile, tmpDcdFile)
      self.o.read_pdb(pdbFile)
      fp = self.o.open_dcd_write(tmpDcdFile)
      self.assertEqual(str(type(fp)),"<type 'SwigPyObject'>")
      os.remove(tmpDcdFile)



   def tearDown(self):
      pass
        
   
   
if __name__ == '__main__': 
   main() 

