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

from sasmol.test_sasmol.util import env

"""
Integration test for sasio.sassubset.Mask.copy_molecule_using_mask

contract:


null test by providing the empty mask

negative test by providing the wrong mask

1-atom/1-frame, mask no atom
1-atom/1-frame, mask 1 atom
1-atom/2-frames, mask frame-0(1 atom)
1-atom/2-frames, mask frame-1(1 atom)

2-aa/1-frames, mask frame-0(no atom)
2-aa/1-frames, mask frame-0(8 atoms)
2-aa/1-frames, mask frame-0(8 atoms)
2-aa/3-frames, mask frame-0(no atom)
2-aa/3-frames, mask frame-0(8 atoms)
2-aa/3-frames, mask frame-0(all atoms)
2-aa/3-frames, mask frame-1(8 atoms)
2-aa/3-frames, mask frame-2(8 atoms)

rna/1-frame, mask frame-0(no atom)
rna/1-frame, mask frame-0(248 atoms)
rna/1-frame, mask frame-0(all atoms)

small protein (crambin)/1-frame, mask frame-0(no atom)
small protein (crambin)/1-frame, mask frame-0(46 atoms)
small protein (crambin)/1-frame, mask frame-0(all atoms)

large protein (groel)/1-frame, mask frame-0(no atom) (Skipped as SASSIE_LARGETEST)
large protein (groel)/1-frame, mask frame-0(7350 atoms) (Skipped as SASSIE_LARGETEST)
large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASSIE_LARGETEST)
large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASSIE_LARGETEST)

a bad pdb, assertRaises
"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset

import numpy

import os

PdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

class Test_sassubset_Mask_copy_molecule_using_mask(MockerTestCase): 
 

   def setUp(self):
      self.o=sasmol.SasMol(0)
      self.o_result=sasmol.SasMol(1)
      self.o_expected=sasmol.SasMol(2)


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


   def generate_expected_pdb(self, o, mask, frame):
      self.o_expected.setAtom([y for (x,y) in zip(mask,self.o.atom()) if x==1])
      self.o_expected.setIndex([y for (x,y) in zip(mask,self.o.index()) if x==1])
      self.o_expected.setName([y for (x,y) in zip(mask,self.o.name()) if x==1])
      self.o_expected.setLoc([y for (x,y) in zip(mask,self.o.loc()) if x==1])
      self.o_expected.setResname([y for (x,y) in zip(mask,self.o.resname()) if x==1])
      self.o_expected.setChain([y for (x,y) in zip(mask,self.o.chain()) if x==1])
      self.o_expected.setResid([y for (x,y) in zip(mask,self.o.resid()) if x==1])
      self.o_expected.setRescode([y for (x,y) in zip(mask,self.o.rescode()) if x==1])
      self.o_expected.setOccupancy([y for (x,y) in zip(mask,self.o.occupancy()) if x==1])
      self.o_expected.setBeta([y for (x,y) in zip(mask,self.o.beta()) if x==1])
      self.o_expected.setSegname([y for (x,y) in zip(mask,self.o.segname()) if x==1])
      self.o_expected.setElement([y for (x,y) in zip(mask,self.o.element()) if x==1])
      self.o_expected.setCharge([y for (x,y) in zip(mask,self.o.charge()) if x==1])
      self.o_expected.setCoor([[y for (x,y) in zip(mask,self.o.coor()[frame]) if x==1]])
      self.o_expected.setNatoms(len(self.o_expected.index()))


   def assert_pdb(self, o1, o2):
      try:
         self.assertEqual(o1.atom(), o2.atom())
         self.assert_list_almost_equal(o1.index(), o2.index())
         self.assertEqual(o1.name(), o2.name())
         self.assertEqual(o1.loc(), o2.loc())
         self.assertEqual(o1.resname(), o2.resname())
         self.assertEqual(o1.chain(), o2.chain())
         self.assert_list_almost_equal(o1.resid(), o2.resid())
         self.assertEqual(o1.rescode(), o2.rescode())
         self.assert_list_almost_equal(o1.coor()[0], o2.coor()[0])
         self.assertEqual(o1.occupancy(), o2.occupancy())
         self.assertEqual(o1.beta(), o2.beta())
         self.assertEqual(o1.segname(), o2.segname())
         self.assertEqual(o1.element(), o2.element())
         self.assertEqual(o1.charge(), o2.charge())
         #self.assertEqual(o1.moltype(), o2.moltype())
      except:
         raise Exception 


   def test_null(self):
      '''
      null
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False
      mask = [0]
      #
      with self.assertRaises(Exception):
         error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)


   def test_wrong(self):
      '''
      negative test by providing the wrong mask
      '''
      #
      try:
         self.o.read_pdb(PdbDataPath+'XXX.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False
      mask = [1,2]
      #
      with self.assertRaises(Exception):
         error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)


   def test_1ATM_mask_none(self):
      '''
	   1-atom/1-frame, mask no atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      frame = 0
      #
      expecting_error = False
      mask = [0]
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1ATM_mask_one(self):
      '''
	   1-atom/1-frame, mask 1 atom
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM.pdb')
      #
      frame = 0
      #
      expecting_error = False 
      mask=[1]
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1ATM_2frames_mask_frame_0(self):
      '''
	   1-atom/2-frames, mask frame-0(1 atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM-1to2.pdb')
      #
      frame = 0
      #
      expecting_error = False 
      mask=[1]

      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1ATM_2frames_mask_frame1(self):
      '''
	   1-atom/2-frames, mask frame-1(1 atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1ATM-1to2.pdb')
      #
      frame = 1
      #
      expecting_error = False 
      mask=[1]
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_1frame_mask_frame0_atomNone(self):
      '''
	   2-aa/1-frames, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B" or resid[i]==12'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_1frame_mask_frame0_8atoms(self):
      '''
	   2-aa/1-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N" or name[i]=="CA" or name[i]=="C" or name[i]=="O"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_1frame_mask_frame0_atoms_all(self):
      '''
	   2-aa/1-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_3frames_mask_frame0_atoms_none(self):
      '''
	   2-aa/3-frames, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_3frames_mask_frame0_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-0(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_3frames_mask_frame0_allAtoms(self):
      '''
	   2-aa/3-frames, mask frame-0(all atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #G
      frame = 0
      basis_filter = 'resid[i]!=5'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_3frames_mask_frame1_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-1(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 1
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_2AAD_3frames_mask_frame2_8atoms(self):
      '''
	   2-aa/3-frames, mask frame-2(8 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'2AAD-1to3.pdb')
      #
      frame = 2
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_rna_1(self):
      '''
      rna/1-frame, mask frame-0(no atom)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==515'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_rna_2(self):
      '''
      rna/1-frame, mask frame-0(248 atoms)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]==20'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_rna_3(self):
      '''
      rna/1-frame, mask frame-0(all atoms)
      '''
      #
      self.o.read_pdb(PdbDataPath+'rna.pdb')
      #
      frame = 0
      basis_filter = 'resid[i]>0'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1CRN_1(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(no atom)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)
      

   def test_1CRN_2(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(46 atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1CRN_3(self):
      '''
      small protein (crambin)/1-frame, mask frame-0(all atoms)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1CRN.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_1(self):
      '''
      large protein (groel)/1-frame, mask frame-0(no atom) (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_2(self):
      '''
      large protein (groel)/1-frame, mask frame-0(7350 atoms) (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]=="N"'
      #
      expecting_error = False
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      print self.o_result.natoms()
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)

   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_1KP8_3(self):
      '''
      large protein (groel)/1-frame, mask frame-0(all atoms) (Skipped as SASSIE_LARGETEST)
	   '''
      #
      self.o.read_pdb(PdbDataPath+'1KP8.pdb')
      #
      frame = 0
      basis_filter = 'name[i]!="B"'
      #
      expecting_error = False 
      errortmp, mask = self.o.get_subset_mask(basis_filter)      
      self.generate_expected_pdb(self.o, mask, frame)
      #
      error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)
      #
      self.assertEqual(len(error)>0, expecting_error)
      self.assert_pdb(self.o_result, self.o_expected)


   def test_1PSI(self):
      '''
      a bad pdb, assertRaises
      '''
      try:
         self.o.read_pdb(PdbDataPath+'1PSI.pdb')
      except Exception:
         pass
      #
      frame = 0
      #
      expecting_error = False 
      mask = [0]*3718   
      #
      with self.assertRaises(Exception):
         error = self.o.copy_molecule_using_mask(self.o_result, mask, frame)



   def tearDown(self):
      pass


if __name__ == '__main__': 
   main() 
