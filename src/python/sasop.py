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
import sasmath

#	SASOP
#
#	12/13/2009	--	initial coding			:	jc
#
#	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Sasop contains the classes and methods to perform the basic
	translation, rotation, and alignment operations on instances of objects 
	described in the sasmol module.

	The methods in this class move entire objects, not pieces
	of single objects (i.e. no intra-object movements)
'''

class Move():

	'''
	This class moves an object.

	masscheck makes sure the mass and center of mass (COM) is current

	translate moves the object to a point in space.
	
	center moves the COM of the object to (0,0,0)

	moveto moves the COM to a point in space (x,y,z)

	'''	

	def masscheck(self,frame):
		if(self._totalmass <=0.0):
			self.calcmass()
		return
	
	def translate(self,frame,value):

		'''
		Simple movement.  It accepts an array of three
		numbers and adds this array to each element.
		It ends by updating the center of mass.   
		'''
		
		self._coor[frame,:,0]=self._coor[frame,:,0]+value[0]
		self._coor[frame,:,1]=self._coor[frame,:,1]+value[1]
		self._coor[frame,:,2]=self._coor[frame,:,2]+value[2]

		self.masscheck(frame)
		self.calccom(frame)

		return

	def center(self,frame):

		'''
		Simple movement.  It moves the center of mass
		to (0.0,0.0,0.0).  The method checks that
		the COM has been calculated first.  
		It ends by updating the center of mass
		'''
	
		self.masscheck(frame)
		self.calccom(frame)
			
		self._coor[frame,:,0]=self._coor[frame,:,0]-self._com[0]
		self._coor[frame,:,1]=self._coor[frame,:,1]-self._com[1]
		self._coor[frame,:,2]=self._coor[frame,:,2]-self._com[2]
		
		self.calccom(frame)
	
		return

	def moveto(self,frame,value):
		'''
		Simple movement.  It moves the center of mass
		to the destination value=[x,y,z].  The method 
		checks that the COM has been calculated first.  
		It ends by updating the center of mass
		'''
		
		self.masscheck(frame)
		self.calccom(frame)
	
		self._coor[frame,:,0]=self._coor[frame,:,0]-self._com[0]+value[0]
		self._coor[frame,:,1]=self._coor[frame,:,1]-self._com[1]+value[1]
		self._coor[frame,:,2]=self._coor[frame,:,2]-self._com[2]+value[2]
		
		self.calccom(frame)
	
		return

	def align(self,frame,coor_sub_2,com_sub_2,coor_sub_1,com_sub_1):

		'''
		Alignment of one object on top of another
		"self" is aligned onto "other" using the basis
		of molecule 2 to align onto the basis of molecule 1
		and the transformation is then done to all the atoms of
		molecule 2

		'''
		self.masscheck(frame)
		self.calccom(frame)
	
		u = sasmath.find_u(coor_sub_1, coor_sub_2)

		tao = numpy.transpose(self.coor()[frame] - com_sub_2)

		error,nat2 = sasmath.matrix_multiply(u,tao)

		ncoor = numpy.transpose(nat2) + com_sub_1

		self._coor[frame,:] = ncoor
 
		return
 
	def rotate(self,frame,axis,theta):
			
		'''
		Simple rotation about the x, y, or z axis.
		
		Note that calcuations are in radians

		'''

		cs=numpy.cos(theta) ; si=numpy.sin(theta)
		if(axis=='x'):
			mat=numpy.array([[1.0,0.0,0.0],[0.0,cs,-si],[0.0,si,cs]])
		elif(axis=='y'):
			mat=numpy.array([[cs,0.0,si],[0.0,1.0,0.0],[-si,0.0,cs]])
		elif(axis=='z'):
			mat=numpy.array([[cs,-si,0.0],[si,cs,0.0],[0.0,0.0,1.0]])

		coordt=self._coor[frame,:].T	
		error,matrix_product = sasmath.matrix_multiply(mat,coordt)
		ncoord = matrix_product.T	
		self._coor[frame,:]=ncoord 
	
		return

	def general_axis_rotate(self,frame,theta,ux,uy,uz):
		'''
        The general rotation of a molecule along an arbitrarily
		given unit axis (ux,uy,uz) by an angle theta.

        Note that calcuations are in radians
        '''

		c11 = numpy.cos(theta)+pow(ux,2)*(1-numpy.cos(theta))
		c12 = ux*uy*(1-numpy.cos(theta))-uz*numpy.sin(theta)
		c13 = ux*uz*(1-numpy.cos(theta))+uy*numpy.sin(theta)
		c21 = uy*ux*(1-numpy.cos(theta))+uz*numpy.sin(theta)
		c22 = numpy.cos(theta)+pow(uy,2)*(1-numpy.cos(theta))
		c23 = uy*uz*(1-numpy.cos(theta))-ux*numpy.sin(theta)
		c31 = uz*ux*(1-numpy.cos(theta))-uy*numpy.sin(theta)
		c32 = uz*uy*(1-numpy.cos(theta))+ux*numpy.sin(theta)
		c33 = numpy.cos(theta)+pow(uz,2)*(1-numpy.cos(theta))

		C = numpy.matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])

		coor = numpy.array(self.coor()[frame]*C)

		self.coor()[frame,:] = coor

		return

	def euler_rotate(self,frame,phi,theta,psi):
		'''
        Rotate the molecule by a euler angle set (phi,theta,psi)

        Note that calcuations are in radians

        '''

		c11 = numpy.cos(theta)*numpy.cos(psi)
		c12 = numpy.cos(phi)*numpy.sin(psi) + numpy.sin(phi)*numpy.sin(theta)*numpy.cos(psi)
		c13 = numpy.sin(phi)*numpy.sin(psi) - numpy.cos(phi)*numpy.sin(theta)*numpy.cos(psi)
		c21 = -numpy.cos(theta)*numpy.sin(psi)
		c22 = numpy.cos(phi)*numpy.cos(psi) - numpy.sin(phi)*numpy.sin(theta)*numpy.sin(psi)
		c23 = numpy.sin(phi)*numpy.cos(psi) + numpy.cos(phi)*numpy.sin(theta)*numpy.sin(psi)
		c31 = numpy.sin(theta)
		c32 = -numpy.sin(phi)*numpy.cos(theta)
		c33 = numpy.cos(phi)*numpy.cos(theta)

		C = numpy.matrix([[c11,c12,c13],[c21,c22,c23],[c31,c32,c33]])
                
		coor = numpy.array(self.coor()[frame]*C)
               
		self.coor()[frame,:] = coor 
                
		return

