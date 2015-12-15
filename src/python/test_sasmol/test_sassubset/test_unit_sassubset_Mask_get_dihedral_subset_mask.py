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
contract:
one residue, mask none
one residue, mask it
three residues, mask none
three residues, mask the first one
three residues, mask the second one
three residues, mask the third one
three residues, mask the first and second ones
three residues, mask the second and third ones
three residues, mask the third and first ones
three residues, mask all residues
500 residues, mask none
500 residues, mask the 300th residue
500 residues, mask the [123,12,90,399,1,89,221,78,91,129]'th residues
500 residues, mask all residues
rna with 1 residue, mask none
rna with 1 residue, mask all
rna with 5 residue, mask none
rna with 5 residue, mask 2,3
rna with 5 residue, mask all
rna with 50 residue, mask all
"""

from unittest import main,skipIf 
from mocker import Mocker, MockerTestCase, ARGS

import sasmol.sasmol as sasmol
import sasmol.sassubset as sassubset

import numpy

import os

class Test_unit_sassubset_Mask_get_dihedral_subset_mask(MockerTestCase): 

   def mock_up(self, Cls_ptch, mthd, mocker, result=None, mmin=0, mmax=None):
      methodToCall = getattr(Cls_ptch,mthd)
      methodToCall(ARGS)
      mocker.result(result)
      mocker.count(mmin, mmax)

   def mock_up_get_dihedral_subset_mask(self, Cls, mocker, natoms, name, resid):
      Cls_ptch = mocker.patch(Cls)
      self.mock_up(Cls_ptch, 'natoms', mocker, natoms)
      self.mock_up(Cls_ptch, 'name', mocker, name)
      self.mock_up(Cls_ptch, 'resid', mocker, resid)
      mocker.replay()

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
 

   def setUp(self):
      self.m=Mocker()
      self.o=sasmol.SasMol(0)


   def test_single_residue_nomask(self):
      '''
      test for a single residue
      nomask for that residue
      '''
      natoms=5
      name=['N','CA','C','O','CB']
      resid=[1]*natoms
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[]
      #
      expected_farray = []
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_single_residue_mask(self):
      '''
      test for a single residue
      mask that residue
      '''
      natoms=5
      name=['N','CA','C','O','CB']
      resid=[1]*natoms
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[1]
      #
      expected_farray = [[1,1,1,0,0]]
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_nomask(self):
      '''
      test for three residue
      nomask for all residues
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[]
      #
      expected_farray = []
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_first_one(self):
      '''
      test for three residues
      mask the first residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[1]
      #
      expected_farray = [[1,1,1,0,0,0]+[1,0,0,0,0,0]+[0]*6]
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_second_one(self):
      '''
      test for three residues
      mask the second residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[2]
      #
      expected_farray = [[0,0,1,0,0,0]+[1,1,1,0,0,0]+[1,0,0,0,0,0]]
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_third_one(self):
      '''
      test for three residues
      mask the second residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[3]
      #
      expected_farray = [[0,0,0,0,0,0]+[0,0,1,0,0,0]+[1,1,1,0,0,0]]
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_first_second(self):
      '''
      test for three residues
      mask the first and second residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[1,2]
      #
      expected_farray = [[1,1,1,0,0,0]+[1,0,0,0,0,0]+[0]*6, [0,0,1,0,0,0]+[1,1,1,0,0,0]+[1,0,0,0,0,0]]

      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_second_third(self):
      '''
      test for three residues
      mask the first and second residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[2,3]
      #
      expected_farray = [[0,0,1,0,0,0]+[1,1,1,0,0,0]+[1,0,0,0,0,0], [0,0,0,0,0,0]+[0,0,1,0,0,0]+[1,1,1,0,0,0]]

      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_third_first(self):
      '''
      test for three residues
      mask the third and first residue
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[3,1]
      #
      expected_farray = [[0,0,0,0,0,0]+[0,0,1,0,0,0]+[1,1,1,0,0,0], [1,1,1,0,0,0]+[1,0,0,0,0,0]+[0]*6]

      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_three_residues_mask_all(self):
      '''
      test for three residues
      mask the all residues
      '''
      natoms=18
      name=['N','CA','C','O','CB','CG']*3
      resid=[1]*6+[2]*6+[3]*6
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[1,2,3]
      #
      expected_farray = [[1,1,1,0,0,0]+[1,0,0,0,0,0]+[0]*6, [0,0,1,0,0,0]+[1,1,1,0,0,0]+[1,0,0,0,0,0], [0,0,0,0,0,0]+[0,0,1,0,0,0]+[1,1,1,0,0,0], ]

      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_500_residues_nomask(self):
      '''
      test for 500 residues
      mask no residue
      '''
      name=['N','CA','C','O','CB','CG1','CG2']
      natom=len(name)
      nres=500
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[]
      #
      expected_farray = []
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_500_residues_mask_number100(self):
      '''
      test for 500 residues
      mask the 300th residue
      '''
      name=['N','CA','C','O','CB','CG1','CG2']
      natom=len(name)
      nres=500
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      nf=300
      flexible_residues=[nf]
      #
      expected_farray = [([0]*natom)*(nf-1) + [0,0,1]+[0]*(natom-3) + [1]*3+[0]*(natom-3) + [1]+[0]*(natom-1) + ([0]*natom)*(nres-nf-2)]

      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_500_residues_mask_100to300(self):
      '''
      test for 500 residues
      mask the 300th residue
      '''
      name=['N','CA','C','O','CB','CG1','CG2']
      natom=len(name)
      nres=500
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=range(100,300)
      #
      expected_farray=[]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0,0,1]+[0]*(natom-3) + [1]*3+[0]*(natom-3) + [1]+[0]*(natom-1) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def test_500_residues_mask_random_10_residues(self):
      '''
      test for 500 residues
      mask the 300th residue
      '''
      name=['N','CA','C','O','CB','CG1','CG2']
      natom=len(name)
      nres=500
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=[123,12,90,399,1,89,221,78,91,129]
      #
      expected_farray=[]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0,0,1]+[0]*(natom-3) + [1]*3+[0]*(natom-3) + [1]+[0]*(natom-1) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   @skipIf(os.environ['SASSIE_LARGETEST']=='n',"I am not testing large files")
   def test_500_residues_mask_all_residues(self):
      '''
      test for 500 residues
      mask the 300th residue
      '''
      name=['N','CA','C','O','CB','CG1','CG2']
      natom=len(name)
      nres=500
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=0
      #
      flexible_residues=range(1,nres-1)
      #
      expected_farray=[]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0,0,1]+[0]*(natom-3) + [1]*3+[0]*(natom-3) + [1]+[0]*(natom-1) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)

   def test_1_residues_rna_mask_none(self):
      '''
      test for rna with 1 residue, mask none
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=1
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      expected_farray = []
      flexible_residues=[]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0]*(natom-1)+[1] + [1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [1]+[0]*(natom-1) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)

   def test_1_residues_rna_mask_all(self):
      '''
      test for rna with 1 residue, mask all
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=1
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      flexible_residues=[0]
      #
      expected_farray = numpy.array([[1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1]],numpy.long)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)

   def test_5_residues_rna_mask_all(self):
      '''
      test for rna with 5 residue, mask none
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=5
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      expected_farray = []
      flexible_residues=[]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0]*(natom-1)+[1] + [1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [1,0,0,1]+[0]*(natom-4) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)

   def test_5_residues_rna_mask_random(self):
      '''
      test for rna with 5 residue, mask 2 and 3
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=5
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      expected_farray = []
      flexible_residues=[2,3]
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0]*(natom-1)+[1] + [1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [1,0,0,1]+[0]*(natom-4) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', map(int,farray.tolist()[0]), '\nexpected_mask:\n',expected_farray[0]
      print 'result_mask:\n', map(int,farray.tolist()[1]), '\nexpected_mask:\n',expected_farray[1]
      self.assert_list_almost_equal(farray, expected_farray)

   def test_5_residues_rna_mask_all(self):
      '''
      test for rna with 5 residue, mask all
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=5
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      expected_farray = []
      flexible_residues=range(1,nres-1)
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0]*(natom-1)+[1] + [1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [1,0,0,1]+[0]*(natom-4) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)

   def test_500_residues_rna_mask_all(self):
      '''
      test for rna with 500 residue, mask all
      '''
      name=['P', 'O1P', 'O2P', "O5'", "C5'", "H5'", "H5''", "C4'", "H4'", "O4'", "C1'", "H1'", 'N1', 'C6', 'H6', 'C2', 'O2', 'N3', 'H3', 'C4', 'O4', 'C5', 'H5', "C2'", "H2''", "O2'", "H2'", "C3'", "H3'", "O3'"]
      natom=len(name)
      nres=5
      name=name*nres
      natoms=natom*nres
      resid=[x/natom for x in range(natoms)]
      resid=numpy.array(resid,numpy.long)
      mtype=1
      #
      expected_farray = []
      flexible_residues=range(1,nres-1)
      for nf in flexible_residues:
         tmp_farray = ([0]*natom)*(nf-1) + [0]*(natom-1)+[1] + [1, 0, 0, 1, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1] + [1,0,0,1]+[0]*(natom-4) + ([0]*natom)*(nres-nf-2)
         expected_farray.append(tmp_farray)
      #
      self.mock_up_get_dihedral_subset_mask(self.o, self.m, natoms, name, resid)
      farray = self.o.get_dihedral_subset_mask(flexible_residues,mtype)
      #
      print 'result_mask:\n', list(farray), '\nexpected_mask:\n',expected_farray
      self.assert_list_almost_equal(farray, expected_farray)


   def tearDown(self):
      self.m.restore()
      self.m.verify()


if __name__ == '__main__': 
   main() 

