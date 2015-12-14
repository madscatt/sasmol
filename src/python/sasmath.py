'''
    SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D. 

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
import sys
import numpy
import math
import matrix_math

#	SASMATH
#
#	12/13/2009	--	initial coding			:	jc
#
#	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Sasmath methods to perform basic mathematical operations
'''

def cross_product(a,b):
	cross = [0]*3
	cross[0] = a[1]*b[2]-a[2]*b[1]
	cross[1] = a[2]*b[0]-a[0]*b[2]
	cross[2] = a[0]*b[1]-a[1]*b[0]
	return numpy.array(cross)

def matrix_multiply(a,b):

	error = []
		
	shape_a = a.shape
	shape_b = b.shape
	dim_a1 = a.shape[0]
	dim_a2 = a.shape[1] ; dim_b1 = b.shape[0]

	try:
		dim_b2 = b.shape[1]
	except:
		dim_b2=1
	
	c = numpy.zeros((dim_a1,dim_b2),numpy.float)
	if( dim_a2 != dim_b1 ):
		message = 'incompatible matrices'
		error.append(message)
		return error,c		

	c = matrix_math.matrix_multiply(a,b,dim_a1,dim_a2,dim_b2)

	return error,c

def find_u(x,y):
	'''
	Method to find the U matrix used to align two molecules
	'''
	b=numpy.zeros(9,numpy.float)
	k=0
	for i in range(3):
		for j in range(3):
			rad=0.0
			for n in range(len(x)):
				rad=rad+y[n,i]*x[n,j]
			numpy.put(b,k,rad)
			k=k+1
	r=numpy.reshape(b,(-1,3))
	r=numpy.mat(r)
	rt=r.T		# transpose of r
	rtr=rt*r	# matrix multiply rt * r
	uk,ak=numpy.linalg.eig(rtr)
	ak=ak.T

	idx = uk.argsort()
	idx = idx[::-1]
	uk = uk[idx]
	ak = ak[idx]

	ak[2]=numpy.cross(ak[0],ak[1])
	rak0=numpy.inner(r,ak[0])
	rak1=numpy.inner(r,ak[1])
	rak0.shape=(1,3)
	rak1.shape=(1,3)
		
	if(uk[0] == 0.0):
		urak0=(10**15)*rak0
	else:
		urak0=(1.0/math.sqrt(abs(uk[0])))*rak0

	if(uk[1] == 0.0):
		urak1=(10**15)*rak1
	else:
		urak1=(1.0/math.sqrt(abs(uk[1])))*rak1

	urak2=numpy.cross(urak0,urak1)
	bk=numpy.zeros((3,3),numpy.float)
	bk[0]=urak0
	bk[1]=urak1
	bk[2]=urak2
	lu=numpy.zeros(9,numpy.float)
	lk=0
	for j in range(3):
		for i in range(3):
			rad=0.0
			for k in range(3):
				rad=rad+bk[k,i]*ak[k,j]
			numpy.put(lu,lk,rad)
			lk=lk+1
	u=numpy.reshape(lu,(-1,3))
	return u

def vec_sub(a,b,c):
	a[0]=b[0]-c[0]
	a[1]=b[1]-c[1]
	a[2]=b[2]-c[2]
	return a
	
def vec_scale(a,b,c):
	a[0]=b*c[0]
	a[1]=b*c[1]
	a[2]=b*c[2]
	return a

def signed_angle(a,b,c):
	'''
	This method calcultes the sign of the angle which is used in the calculation of a dihedral angle.
	As such, this method is not validated for other uses.  It will fail if the basis atoms for the
	dihedral (atom 2 and atom 3) overlap.
	'''

	ada = (a[0]*a[0] + a[1]*a[1] + a[2]*a[2])
	bdb = (b[0]*b[0] + b[1]*b[1] + b[2]*b[2])
	
	if(ada*bdb <= 0.0) :
		return 180.0
	else:
		adb = (a[0]*b[0] + a[1]*b[1] + a[2]*b[2])
		try:	
			argument = adb/math.sqrt(ada*bdb)
			angle = (180.0/math.pi) * math.acos(argument)
		#angle = (180.0/math.pi) * math.acos(adb/math.sqrt(ada*bdb))
		except:
			#print '>>> exception averted: arg = ',argument
			print ';',
			return 180.0
	
	cp = cross_product(a,b)
	dp = cp[0]*c[0] + cp[1]*c[1] + cp[2]*c[2]
	sign = cmp(dp, 0.0)

	return sign*angle

def dihedral_angle(a1,a2,a3,a4):
	
	r1=numpy.zeros(3,numpy.float)
	r2=numpy.zeros(3,numpy.float)
	r3=numpy.zeros(3,numpy.float)
	r4=numpy.zeros(3,numpy.float)

#	vec_sub(r1,a1,a2)
#	vec_sub(r2,a3,a2)
#	vec_scale(r3,-1,r2)
#	vec_sub(r4,a4,a3)

	r1 = a1 - a2
	r2 = a3 - a2
	r3 = -1.0*r2
	r4 = a4 - a3

	n1=cross_product(r1,r2)
	n2=cross_product(r3,r4)

	dihedral_angle = signed_angle(n1,n2,r2)

	return dihedral_angle

def calc_angle(coor1,coor2,coor3):
	'''
	Method to calculate the angle between three atoms
	'''
	
	u = coor1-coor2 ; v = coor3-coor2

	norm_u = math.sqrt(sum(u*u)) ; norm_v = math.sqrt(sum(v*v))
	c = numpy.dot(u,v)/(norm_u*norm_v) # -> cosine of the angle
	angle = math.acos(c) 

	return angle

