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
from mocker import Mocker, MockerTestCase, ANY, ARGS, KWARGS
import sasmol.sasmol as sasmol

import numpy

import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_sascalc_Prop_calcrmsd(MockerTestCase): 

    def setUp(self):
        self.o1=sasmol.SasMol(0)
        self.o2=sasmol.SasMol(0)

    def calc_exp(self):
        c1 = numpy.array((self.o1.coor()[0]),floattype)
        c2 = numpy.array((self.o2.coor()[0]),floattype)
        return numpy.sqrt(numpy.sum((c1-c2)**2)/len(c1))

    def test_null(self):
        self.o1.setCoor(numpy.zeros((1,0,3),floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.zeros((1,0,3),floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        self.assertTrue(numpy.isnan(result_rmsd))

    def test_one_overlap_atom(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = 0.0
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_one_nonoverlap_atom(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[4.0, 5.0, 6.0]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = 3.0*numpy.sqrt(3.0)
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_two_atoms(self):
        self.o1.setCoor(numpy.array([[[7.0, 8.0, 9.0],[1.0, 3.0, 5.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[12.0, 53.0, 67.0],[76.0, 87.0, 96.0]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = self.calc_exp()
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_six_atoms(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, 53.0, 67.0],[76.0, 87.0, 96.0],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = self.calc_exp()
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_six_atoms_inf1(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, util.HUGE, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, util.HUGE, 67.0],[76.0, 87.0, 96.0],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = util.INF
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_6_atoms_inf2(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.INF, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, 53.0, 67.0],[76.0, 87.0, util.INF],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = util.INF
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_6_atoms_nan(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.NAN, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, 53.0, 67.0],[76.0, 87.0, util.NAN],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        self.assertTrue(numpy.isnan(result_rmsd))

    def test_6_atoms_tiny(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.TINY, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, 53.0, 67.0],[76.0, 87.0, util.TINY],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = self.calc_exp()
        self.assertAlmostEqual(expected_rmsd, result_rmsd)

    def test_6_atoms_zero(self):
        self.o1.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.ZERO, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o1.setNatoms(len(self.o1._coor[0]))
        self.o2.setCoor(numpy.array([[[2.0, 12.0, 35.0],[12.0, 53.0, 67.0],[76.0, 87.0, util.ZERO],[12.0, 33.0, 52.0],[2.3, 4.3, 6.8],[0.0,  22.5,33.6]]],floattype))
        self.o2.setNatoms(len(self.o2._coor[0]))
        result_rmsd  = self.o1.calcrmsd(self.o2)
        expected_rmsd = self.calc_exp()
        self.assertAlmostEqual(expected_rmsd, result_rmsd)


    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

