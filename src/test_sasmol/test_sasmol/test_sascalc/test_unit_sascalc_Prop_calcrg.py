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

from sassie.core_testing.util import env, util

from unittest import main 
from mocker import Mocker, MockerTestCase
from sassie.sasmol import sasmol

import numpy


import warnings; warnings.filterwarnings('ignore')

import os
floattype=os.environ['SASSIE_FLOATTYPE']

class Test_sascalc_Prop_calcrg(MockerTestCase): 

    def setUp(self):
        self.o=sasmol.SasMol(0)
        self.tol = 3

    def calc_exp(self):
        self.o.calccom(0)
        coor = numpy.array((self.o.coor()[0]),floattype)
        return numpy.sqrt(numpy.sum((coor-self.o.com())**2)/len(coor))

    def test_one_null(self):
        with self.assertRaises(Exception):
          result_rg  = self.o.calcrg(0)

    def test_one_atom(self):
        self.o.setElement(['C'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o.element()))
        result_rg  = self.o.calcrg(0)
        expected_rg = 0.0
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_two_atoms(self):
        self.o.setElement(['C', 'AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o.element()))
        result_rg  = self.o.calcrg(0)
        expected_rg = 3.3265416
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_six_atoms_duplicate(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calcrg(0)
        expected_rg = self.calc_exp()
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_six_atoms_duplicate_inf1(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.HUGE, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calcrg(0)
        expected_rg = util.INF
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_six_atoms_duplicate_inf2(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.INF, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calcrg(0)
        self.assertTrue(numpy.isnan(result_rg))

    def test_six_atoms_duplicate_nan(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.NAN, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calcrg(0)
        self.assertTrue(numpy.isnan(result_rg))

    def test_six_atoms_duplicate_tiny(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.TINY, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))
        result_rg  = self.o.calcrg(0)
        expected_rg = self.calc_exp()
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_six_atoms_duplicate_zero(self):
        self.o.setElement(['C','O','C','1H','KR','AG'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[util.ZERO, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))        
        result_rg  = self.o.calcrg(0)
        expected_rg = self.calc_exp()
        self.assertAlmostEqual(expected_rg, result_rg, self.tol)

    def test_wrong(self):
        self.o.setElement(['X','M'])
        self.o.setCoor(numpy.array([[[1.0, 2.0, 3.0],[7.0, 8.0, 9.0]]],floattype))
        self.o.setTotalmass(0.0)
        self.o.setNatoms(len(self.o._element))        
        result_rg  = self.o.calcrg(0)
        import math
        self.assertTrue(math.isnan(result_rg))



    def tearDown(self):
        pass

if __name__ == '__main__': 
   main() 

