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
import os
import sys
import string
import locale
import struct
import numpy
import time
import sasview_vmd

#	SASVIEW
#
#	11/27/2013	--	initial coding			:	jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Sasview is the main module that contains the base classes that 
	read and write atomic information from and to external viewers
	such as vmd and others in the future (PyMol, etc.).
	
	To view the coordinates the viewer needs to be open and controlled
	external to sasview / sasmol.  These methods merely allow 
	coordinates to be passed between programs.

	See the following sites for the IMD format:

	http://www.ks.uiuc.edu/Research/vmd/
	http://www.ks.uiuc.edu/Research/namd/

	These classes are accessed by the SasAtm class found in
	the sasmol module.

'''

class View(object):

	def __init__(self,filename,flag):
		pass

	def check_error(self,error):

		if(len(error)>0):
			print error
			sys.exit()

		return


	def send_coordinates_to_vmd(self,port,flag):
		'''
		This method opens a socket to send and receive coordinates
		by calling a pre-compiled C module (sasview_vmd).
		'''
		natoms = self._coor[0,:,0].shape[0]
		frame = 0		
		tx=self._coor[frame,:,0].astype(numpy.float32)	
		ty=self._coor[frame,:,1].astype(numpy.float32)	
		tz=self._coor[frame,:,2].astype(numpy.float32)	

		result = sasview_vmd.send_coordinates_to_vmd(tx,ty,tz,port,flag)

		print 'result = ',result

		return 

