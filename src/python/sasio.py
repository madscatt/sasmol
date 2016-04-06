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
#### NOTE: INDENTATION FIXED : 12/11/2015 JEC

import os
import sys
import string
import locale
import struct
import numpy
import time
import dcdio

#	SASIO
#
#	12/5/2009	--	initial coding			:	jc
#	12/10/2009	--	doc strings 			:	jc
#	01/01/2011	--	added dcdio wrappers		:	jc
#
#LC	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Sasio is the main module that contains the base classes that 
	read and write atomic information from and to the hard disk,
	and (eventually) deal with logging of input parameters and
	project runs.

	The methods in class Files are used to read and write data to
	the Charmm/Xplor binary data format (DCD) and textual protein
	data bank (PDB) format.  

	See the following sites for the DCD format:

	http://www.ks.uiuc.edu/Research/vmd/plugins/molfile/dcdplugin.html
	http://www.lrz-muenchen.de/~heller/ego/manual/node93.html

	The methods (read_pdb, write_pdb) are used to read and write data
	using the PDB file format as described at the following
	web-site:

	http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

	These classes are accessed by the SasAtm class found in
	the sasmol module.

'''

class Files(object):

    def __init__(self,filename,flag):
        pass

    def open_dcd_read(self,filename):
        '''
        This method opens a file to read in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        filepointer=dcdio.open_dcd_read(filename)

        num_fixed=0
        result = 1

        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(filepointer)
        if(readheaderresult!=0):
            print 'failed to read header'
            print 'readheaderresult = ',readheaderresult	
	
        dcdfile = [filepointer,nnatoms,nset,reverseEndian,charmm]
	
        return dcdfile

    def open_dcd_write(self,filename):
        '''
        This method opens a file and writes the headerfile in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).

        This function will OVERWRITE a file with the same name without prompting.
        '''

        filepointer=dcdio.open_dcd_write(filename)

        self.write_dcd_header(filepointer,1)

        return filepointer

    def write_dcd_header(self,filepointer,nset):
        '''
        This method writes the headerfile in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''
        filename=" "
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0
        headerresult=dcdio.write_dcdheader(filepointer,filename,natoms,nset,istart,nsavc,delta)

        return 

    def write_dcd_step(self,filepointer,frame,step):
        '''
        This method writes a single step in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        tx=self._coor[frame,:,0].astype(numpy.float32)	
        ty=self._coor[frame,:,1].astype(numpy.float32)	
        tz=self._coor[frame,:,2].astype(numpy.float32)	

        stepresult=dcdio.write_dcdstep(filepointer,tx,ty,tz,step)

        return

    def write_dcd_frames(self, filename, start, end):
        '''
        This method writes a single step or multiple frames
        in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''

        outfile=dcdio.open_dcd_write(filename)
        nset = end-start
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0

        headerresult=dcdio.write_dcdheader(outfile,filename,natoms,nset,istart,nsavc,delta)

        i = 0
        for frame in xrange(start,end):
            print ".",
            sys.stdout.flush()

            tx=self._coor[frame,:,0].astype(numpy.float32)	
            ty=self._coor[frame,:,1].astype(numpy.float32)	
            tz=self._coor[frame,:,2].astype(numpy.float32)	

            stepresult=dcdio.write_dcdstep(outfile,tx,ty,tz,i+1)
            i += 1

        result = dcdio.close_dcd_write(outfile)

        return

    def close_dcd_write(self,filepointer):
        '''
        This method closes a dcd file.
        '''
        result = dcdio.close_dcd_write(filepointer)
        print 'result = ',result
	
        return
	
    def close_dcd_read(self,filepointer):
        '''
        This method closes a dcd file.
        '''
        result = dcdio.close_dcd_read(filepointer)
        print 'result = ',result
	
        return

    def write_dcd(self,filename):
        '''
        This method writes data in the Charmm/Xplor data format.
        by calling a pre-compiled C module (dcdio).
        '''
		
        outfile=dcdio.open_dcd_write(filename)

        nset = self._coor[:,0,0].shape[0]
        natoms = self._coor[0,:,0].shape[0]
        istart = 0 ; nsavc = 1 ; delta = 1.0
		
        headerresult=dcdio.write_dcdheader(outfile,filename,natoms,nset,istart,nsavc,delta)

        for frame in xrange(nset):
            print ".",
            sys.stdout.flush()

            tx=self._coor[frame,:,0].astype(numpy.float32)	
            ty=self._coor[frame,:,1].astype(numpy.float32)	
            tz=self._coor[frame,:,2].astype(numpy.float32)	

            stepresult=dcdio.write_dcdstep(outfile,tx,ty,tz,frame+1)

        result = dcdio.close_dcd_write(outfile)

        return

    def read_single_dcd_step(self,filename,frame):
        '''
        This method reads a single dcd step in the Charmm/Xplor data format.
	
        The method simply reads all frames up until frame and then assigns
        coordinates to the last frame (no seek option is utilizied)

        '''

        infile=dcdio.open_dcd_read(filename)
        num_fixed=0
        result = 1

        print 'calling read dcd header'
        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(infile)
        if(readheaderresult!=0):
            print 'failed to read header'
            print 'readheaderresult = ',readheaderresult	

        print 'done with read dcd header'

        coor=numpy.zeros((1,nnatoms,3),numpy.float)	
	
        tx=numpy.zeros(nnatoms,dtype=numpy.float32)
        ty=numpy.zeros(nnatoms,dtype=numpy.float32)
        tz=numpy.zeros(nnatoms,dtype=numpy.float32)

        first = 1  # since num_fixed = 0 ; the "first" variable is inconsequential
		
        print 'calling read_dcdstep'
        for i in xrange(frame):	
            result=dcdio.read_dcdstep(infile,tx,ty,tz,num_fixed,first,reverseEndian,charmm)

        print 'back from read_dcdstep'
        print 'result = ',result

        coor[0,:,0]=tx.astype(numpy.float) ; coor[0,:,1]=ty.astype(numpy.float) ; coor[0,:,2]=tz.astype(numpy.float)
		
        result = dcdio.close_dcd_read(infile)
        self._coor=numpy.array(coor)
	
        if(result!=0):
            print 'failed to read coordinates'	
            print 'result = ',result

        return
	
    def read_dcd_step(self,dcdfile,frame,**kwargs):
        '''
        This method reads a single dcd step in the Charmm/Xplor data format.
        '''
        num_fixed=0

        filepointer = dcdfile[0]
        nnatoms = dcdfile[1]
        reverseEndian = dcdfile[3]
        charmm = dcdfile[4]

        tx=numpy.zeros(nnatoms,dtype=numpy.float32)
        ty=numpy.zeros(nnatoms,dtype=numpy.float32)
        tz=numpy.zeros(nnatoms,dtype=numpy.float32)

        result=dcdio.read_dcdstep(filepointer,tx,ty,tz,num_fixed,frame,reverseEndian,charmm)

        self._coor[0,:,0]=tx.astype(numpy.float) ; self._coor[0,:,1]=ty.astype(numpy.float) ; self._coor[0,:,2]=tz.astype(numpy.float)
	
        if len(kwargs) < 1:
            sys.stdout.write('.',)
        elif not kwargs['no_print']: 
            try:
                sys.stdout.write('.',)
            except:
                pass

        return

    def read_dcd(self,filename):

        '''
        This method reads data in the Charmm/Xplor data format.
        '''

        infile=dcdio.open_dcd_read(filename)

        nnatoms=0 ; nset=0 ; istart=0 ; nsavc=0 ; delta=0.0
        namnf=0 ; freeindexes=[] ; reverseEndian=0 ; charmm=0

        readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(infile)
        coor=numpy.zeros((nset,nnatoms,3),numpy.float)	
	
        num_fixed=0 
        result=1

        sum=0.0
        for i in xrange(nset):
            print '.',
            sys.stdout.flush()
            read_start_time=time.time()
            tx=numpy.zeros(nnatoms,dtype=numpy.float32)
            ty=numpy.zeros(nnatoms,dtype=numpy.float32)
            tz=numpy.zeros(nnatoms,dtype=numpy.float32)
		
            result=dcdio.read_dcdstep(infile,tx,ty,tz,num_fixed,i,reverseEndian,charmm)
            read_end_time=time.time()
	
            sum+=read_end_time-read_start_time

            coor[i,:,0]=tx.astype(numpy.float) ; coor[i,:,1]=ty.astype(numpy.float) ; coor[i,:,2]=tz.astype(numpy.float)
	
        result = dcdio.close_dcd_read(infile)
        self._coor=numpy.array(coor)

        print

        return

    def print_error(self,name,my_message):

        error = []
        message = '\nELEMENT NAME NOT IN CHARMM DATABASE AND SINGLE CHARACTER OPTION NOT APPLICABLE\n'
        message += '\n'+name+'\n'
        message +=  my_message+'\n'
        message += ' stopping now: name = '+name+'\n'
        error.append(message)

        return error

    def check_error(self,error):

        if(len(error)>0):
            print error
            sys.exit()

        return

    def element_filter(self):
        '''	
        This function filters the PDB file to sraighten
        out various things (add to this as you go along).
        
        Filter element list... 
        '''
        unique_elements = []
        error = []
        single_atom_names = ['H','F','B','D','C','N','O','S','P','I','K','U','V','W','Y']

        for i in range(len(self._name)):
            name = self._name[i].upper()
            resname = self._resname[i].upper()
            if (self._element[i]=='' or self._element[i]==' ' or self._element[i]=='  '):
                natoms = len(self._name[i])
                error,self._element[i] = self.get_element(name,resname)
                self.check_error(error)
#
###	OPEN	Error exception handling stub
#
            if(self._element[i] not in unique_elements): 
                unique_elements.append(self._element[i])

        return unique_elements

    def get_element(self,name,resname):
        '''
        Get elment from charmm27 atom name list plus extras from periodic table
        and resolve name conflicts	
        '''

        error = []
        element_name = ''
	
        hydrogen,carbon,nitrogen,oxygen,sulfur,phosphorus,other = self.charmm_names()	
	
        conflict_atoms = ['CD','CE','HE','HG','NE','ND','NB','PB','PA']
	
        if(name in conflict_atoms and (name == resname)):
            element_name = name	
        elif(name in hydrogen): element_name = 'H'
        elif(name in carbon): element_name = 'C'
        elif(name in nitrogen): element_name = 'N'
        elif(name in oxygen): element_name = 'O'
        elif(name in sulfur): element_name = 'S'
        elif(name in phosphorus): element_name = 'P'
        elif(name == '1H'): element_name = name
        elif(name == '2H'): element_name = "D"
        elif(name == "SOD"): element_name = 'NA'
        elif(name == "POT"): element_name = 'K'
        elif(name == "CAL"): element_name = 'CA'
        elif(name == "CLA"): element_name = 'CL'
        elif(name == "CES"): element_name = 'CS'
        elif(name == "F2'"): element_name = 'F'
        elif(name in other): element_name = name
        elif(len(name)>0):
            if(name[0].isalpha()):
                element_name = name[0]
            elif(name[0].isdigit() and len(name)>1):
                if(name[1].isalpha()):
                    element_name = name[1]
                else:
                    my_message = '\nFIRST AND SECOND CHARACTER IN ATOM NAME CAN NOT BE A NUMBER\n'
                    error = self.print_error(name,my_message)
            else:
                my_message = '\nFIRST CHARACTER IN ATOM NAME CAN NOT BE A NUMBER\n'
                error = self.print_error(name,my_message)
        else:
            my_message = '\nNO CHARACTERS FOUND FOR ATOM NAME\n'
            error = self.print_error(name,my_message)
	
        return error,element_name

    def write_pdb(self,filename,frame,flag,**kwargs):

        '''
        This method writes the PDB file
        '''

        debug=0
        result=1
        conect = False
	
        if(flag=='w' or flag=='W'):
            infile=open(filename,'w')
        elif(flag=='a' or flag=='A'):				
            infile=open(filename,'a')
	
        if 'conect' in kwargs:
            conect = kwargs['conect']
		
        if 'model' in kwargs:
            this_frame = kwargs['model']
            infile.write("MODEL "+str(this_frame)+"\n")

        if(debug==1):
            infile.write("          1         2         3         4         5         6         7         \n")
            infile.write("01234567890123456789012345678901234567890123456789012345678901234567890123456789\n")

        for i in range(len(self._atom)):
			
            this_index = self._index[i]
	
            this_resid = self._resid[i]

            if(this_index > 99999):
                this_index = '99999'
            elif(this_index < -9999):
                this_index = '-9999'
            else:
                this_index = str(this_index)

            if(this_resid > 9999):
                this_resid = '9999'
            elif(this_resid < -999):
                this_resid = '-999'

            try:
                sx = "{0:8.3f}".format(float(self._coor[frame,i,0]))[:8]
                sy = "{0:8.3f}".format(float(self._coor[frame,i,1]))[:8]
                sz = "{0:8.3f}".format(float(self._coor[frame,i,2]))[:8]

                infile.write("%-6s%5s %-4s%1s%-4s%1s%4s%1s   %8s%8s%8s%6s%6s      %-4s%2s%2s\n" % (self._atom[i],this_index,self._name[i],self._loc[i],self._resname[i],self._chain[i],this_resid,self._rescode[i],sx,sy,sz,self._occupancy[i],self._beta[i],self._segname[i],self._element[i],self._charge[i]))
            except:
                print '\n>>>> ERROR IN WRITE_PDB <<<<\n'
                print '>> i = ',i
                print 'atom = ',self._atom[i],' : type = ',type(self._atom[i])
                print 'index = ',this_index,' : type = ',type(this_index)
                print 'name = ',self._name[i],' : type = ',type(self._name[i])
                print 'loc = ',self._loc[i],' : type = ',type(self._loc[i])
                print 'resname = ',self._resname[i],' : type = ',type(self._resname[i])
                print 'chain = ',self._chain[i],' : type = ',type(self._chain[i])
                print 'resid = ',self._resid[i],' : type = ',type(self._resid[i])
                print 'rescode = ',self._rescode[i],' : type = ',type(self._rescode[i])
                print 'coor_x = ',self._coor[frame,i,0],' : type = ',type(self._coor[frame,i,0])
                print 'coor_y = ',self._coor[frame,i,1],' : type = ',type(self._coor[frame,i,1])
                print 'coor_z = ',self._coor[frame,i,2],' : type = ',type(self._coor[frame,i,2])
                print 'occupancy = ',self._occupancy[i],' : type = ',type(self._occupancy[i])
                print 'beta = ',self._beta[i],' : type = ',type(self._beta[i])
                print 'segname = ',self._segname[i],' : type = ',type(self._segname[i])
                print 'element = ',self._element[i],' : type = ',type(self._element[i])
                print 'charge = ',self._charge[i],' : type = ',type(self._charge[i])

        # TODO: Check with Joseph to see if logic acceptable -
        # i.e. is 'final' always accompanied by 'model'?
        if ('final' in kwargs) or ('model' not in kwargs):

            if conect:
                conect_lines = self.create_conect_pdb_lines()
                for line in conect_lines:
                    infile.write(line + '\n')

            infile.write("END\n")

        else:

            infile.write("ENDMDL\n")

        # if 'model' in kwargs:
        #     infile.write("ENDMDL\n")
        # else:
        #     if conect:
        #         conect_lines = self.create_conect_pdb_lines()
        #         for line in conect_lines:
        #             infile.write(line + '\n')
        #
        #     infile.write("END\n")
        #
        # if 'final' in kwargs:
        #     infile.write("END\n")

        infile.close()
    
        return result

    def get_resnames(self):
        '''
        This method holds names of residues to use to set moltype.  Based on Charmm 27 naming.
        '''

        protein_resnames=['ALA','ARG','ASN','ASP','CYS','GLN','GLU','GLY','HIS','HSD','HSE','HSP','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']
	
        dna_resnames=['NUSA','NUSG','NUSC','NUSU','DA','DG','DC','DT','ADE','GUA','CYT','THY']

        rna_resnames=['RNUS','RNUA','RUUG','RNUC','A', 'C', 'G', 'U','ADE','CYT','GUA','URA']

        nucleic_resnames = ['GUA','ADE','CYT','THY','URA','G', 'A', 'C', 'T', 'U','DA','DG','DC','DT']

        water_resnames=['TIP3','SPCE','TIP','SPC','TIP4','TP3M'] 
			
        return protein_resnames,dna_resnames,rna_resnames,nucleic_resnames,water_resnames

    def initialize_children(self):

        '''
            This initializes the following "children" with their
            masks already defined to the "parent" molecule

            names() 		: names_mask()
            resnames()		: resnames_mask()
            resids()		: resids_mask()
            chains()		: chains_mask()
            segnames()		: segnames_mask()
            occupancies()	: occupancies_mask()
            betas()		: betas_mask()
            elements()		: elements_mask()
			
            The objects on the left contain the unique values and
            the objects on the right contain the masks that have
            the indices to extract the information for each unique
            value from the parent molecule.

        '''

        number_of_names = self.number_of_names() ; unique_names = self.names()
        number_of_resnames = self.number_of_resnames() ; unique_resnames = self.resnames()
        number_of_resids = self.number_of_resids() ; unique_resids = self.resids()
        number_of_chains = self.number_of_chains() ; unique_chains = self.chains()
        number_of_occupancies = self.number_of_occupancies() ; unique_occupancies = self.occupancies()
        number_of_betas = self.number_of_betas() ; unique_betas = self.betas()
        number_of_elements = self.number_of_elements() ; unique_elements = self.elements()
        number_of_segnames = self.number_of_segnames() ; unique_segnames = self.segnames()
        natoms = self.natoms()
        name = self.name() ; resname = self.resname() ; resid = self.resid() ; chain = self.chain()
        occupancy = self.occupancy() ; beta = self.beta() ; element = self.element() ; segname = self.segname()

        names = [] ; resnames = [] ; resids = [] ; chains = [] 
        occupancies = [] ; betas = [] ; elements = [] ; segments = []

        names_mask = numpy.zeros((number_of_names,natoms),numpy.int)
        resnames_mask = numpy.zeros((number_of_resnames,natoms),numpy.int)
        resids_mask = numpy.zeros((number_of_resids,natoms),numpy.int)
        chains_mask = numpy.zeros((number_of_chains,natoms),numpy.int)
        occupancies_mask = numpy.zeros((number_of_occupancies,natoms),numpy.int)
        betas_mask = numpy.zeros((number_of_betas,natoms),numpy.int)
        elements_mask = numpy.zeros((number_of_elements,natoms),numpy.int)
        segnames_mask = numpy.zeros((number_of_segnames,natoms),numpy.int)

        nresid = 0

        for i in xrange(natoms):

            names_mask[unique_names.index(name[i])][i] = 1
            resnames_mask[unique_resnames.index(resname[i])][i] = 1
            resids_mask[unique_resids.index(resid[i])][i] = 1
            chains_mask[unique_chains.index(chain[i])][i] = 1
            occupancies_mask[unique_occupancies.index(occupancy[i])][i] = 1
            betas_mask[unique_betas.index(beta[i])][i] = 1
            elements_mask[unique_elements.index(element[i])][i] = 1
            segnames_mask[unique_segnames.index(segname[i])][i] = 1

        self._names_mask = names_mask
        self._resnames_mask = resnames_mask
        self._resids_mask = resids_mask
        self._chains_mask = chains_mask
        self._occupancies_mask = occupancies_mask
        self._betas_mask = betas_mask
        self._elements_mask = elements_mask
        self._segnames_mask = segnames_mask

        return


    def check_for_all_zero_columns(self, coor, frame=0):
        '''
            Make sure there are no two all zero columns in coordinates
            This is needed for the alignment module 
        '''
        SMALL = 1.0E-10
        x0,y0,z0=False,False,False
        if sum(coor[frame][:,0]**2)<SMALL:
            x0=True
        if sum(coor[frame][:,1]**2)<SMALL:
            y0=True
        if sum(coor[frame][:,2]**2)<SMALL:
            z0=True
        if x0:
            coor[frame][0][0]=SMALL
        if y0:
            coor[frame][0][1]=SMALL
        if z0:
            coor[frame][0][2]=SMALL

        return

    def read_pdb(self,filename,**kwargs):
        '''
        This method reads a PDB file.  
        '''
        debug=0
        result=1

        protein_resnames,dna_resnames,rna_resnames,nucleic_resnames,water_resnames = self.get_resnames()

        ### LONG	Need to properly import header information
	
        # see: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html

        fastread = False
        pdbscan = False
        printme = True

        if 'verbose' in kwargs:
            printme = kwargs['verbose']
			
        if 'fastread' in kwargs:
            fastread = kwargs['fastread']
            if(fastread):
                printme = False

        if 'pdbscan' in kwargs:
            pdbscan = kwargs['pdbscan']
	
        infile=open(filename,'r').readlines()

        if(printme): print 'reading filename: ',filename
		
        atom=[] ; index=[] ; original_index=[] ; name=[] ; loc=[] ; resname=[] ; chain=[] ; resid=[] ; rescode=[]
        x=[] ; y=[] ; z=[] 
        occupancy=[] ; beta=[] ; segname=[] ; element=[] ; charge=[] ; moltype=[] ; conect = {}
        residue_flag = [] ; original_resid=[] ; header = []
	
        # first check to see how many frames are in the file

        num_model = 0
        num_endmdl = 0
        num_end = 0
        
        num_frames = 1
        
        count_index = 0

        num_counts_this_model = 0 # number of atoms encompassed by "MODEL" and "ENDMDL" lines
        num_counts_this_end = 0 # number of atoms before the first "END" line or between two consecutive "END" lines
        num_counts_per_model = [] # array of number of atoms encompassed by "MODEL" and "ENDMDL" lines
        num_counts_per_end = [] # array of number of atoms before the first "END" line or between two consecutive "END" lines
        modelON = False # Flag to indicate that a "MODEL" frame is being read
        for i in range(len(infile)):
            lin=infile[i]
            lins=string.split(infile[i])
            if (len(lins)==0):
                if ((i+1)==len(infile)):
                    lins=['']
                    if modelON:
                        raise Exception, 'There should be an ENDMDL pairing with MODEL'
                    else:
                        continue

            record_name = string.strip(lin[0:6])
#
### 	OPEN 	need to re-factor the exception statements to a uniform reporting mechanism
#
            try:
                if(lins[0]=='MODEL'):
                    if modelON:
                        raise Exception, 'Encountered two consecutive MODEL lines' 
                    if (num_counts_this_model != 0):
                        raise Exception, 'There should not be atoms after ENDMDL and before MODEL lines'
                    modelON = True
                elif(lins[0]=='ENDMDL'):
                    if not modelON:
                        raise Exception, 'Encountered two consecutive ENDMDL lines'
                    modelON = False
                    num_counts_per_model.append(num_counts_this_model)
                    num_counts_this_model = 0
                elif(lins[0]=='END'):
                    num_counts_per_end.append(num_counts_this_end)
                    num_counts_this_end = 0
            except:
                pass
			#
            if((record_name == 'ATOM' or record_name == 'HETATM')):
                count_index += 1
                num_counts_this_model += 1
                num_counts_this_end += 1
		#
		#

        if ( (len(num_counts_per_end)==0) and (len(num_counts_per_model)!=0) ):
            raise Exception, 'According to Protein Data Bank Contents Guide, END line must appear in each coor entry'
        if (len(num_counts_per_model)!=0 and (len(num_counts_per_end)>1 or sum(num_counts_per_model)!=sum(num_counts_per_end))):
            if(printme): print num_counts_per_model,num_counts_per_end
            raise Exception, 'Only one terminating END line is allowed for pdb entries with multiple MODEL'
		#
        if (len(num_counts_per_model)>0):
            num_frames = len(num_counts_per_model)
            num_atoms = num_counts_per_model[0]
            if not all(x == num_atoms for x in num_counts_per_model):
                raise Exception, 'number of atoms per frame is not equal'
        elif (len(num_counts_per_end)>0):
            num_frames = len(num_counts_per_end)
            num_atoms = num_counts_per_end[0]
            if not all(x == num_atoms for x in num_counts_per_end):
                raise Exception, 'number of atoms per frame is not equal'
        elif ( (len(num_counts_per_model)==0) and (len(num_counts_per_end)==0) ):
            num_frames = 1
            num_atoms = num_counts_this_model
        else:
            raise Exception, 'unexpected error!'

        if(printme): print 'num_atoms = ',num_atoms

        if(printme): print '>>> found ',num_frames,' model(s) or frame(s)'			
        this_frame = 1		

        true_index = 0
	
        coor=numpy.zeros((num_frames,num_atoms,3),numpy.float)	

        unique_names = [] ; unique_resnames = [] ; unique_resids = [] ; unique_chains = [] 
        unique_occupancies = [] ; unique_betas = [] ; unique_segnames = [] ; unique_moltypes = []

        for i in range(len(infile)):
            lin=infile[i]
            lins=string.split(infile[i])

            record_name = string.strip(lin[0:6])

            if((record_name == 'ATOM' or record_name == 'HETATM') and this_frame == 1):
                true_index += 1
                atom.append(string.strip(lin[0:6]))		#	1-6		record name	
                original_index.append(lin[6:11])				#	7-11		atom serial number
                index.append(str(true_index))   	        #   set index so that > 99,999 atoms can be read and counted
                this_name = string.strip(lin[12:16])		#	13-16		atom name
                name.append(string.strip(lin[12:16]))		#	13-16		atom name
                if pdbscan:
                    loc.append(lin[16])
                else:
                    loc.append(' ')
                this_resname = string.strip(lin[17:21])	#	18-20		residue name
                resname.append(string.strip(lin[17:21]))	#	18-20		residue name
                this_chain = lin[21]				#	22		chain identifier
                chain.append(lin[21])				#	22		chain identifier
                this_resid = locale.atoi(lin[22:26])			#	23-26		residue sequence number
                original_resid.append(lin[22:26])			#	23-26		residue sequence number
                resid.append(lin[22:26])			#	23-26		residue sequence number
                rescode.append(lin[26])				#	27		code for insertion of residues
                x.append(lin[30:38])				#	31-38		Real(8.3) X: angstroms	
                y.append(lin[38:46])				#	39-46		Real(8.3) Y: angstroms	
                z.append(lin[46:54])				#	47-54		Real(8.3) Z: angstroms	
				
                residue_flag.append(False)

                if not pdbscan:

                    try:	
                        occupancy.append(string.strip(lin[54:60]))      #	55-60		occupancy
                        this_occupancy =string.strip(lin[54:60])      #	55-60		occupancy
                        if(occupancy[-1] == ''):
                            occupancy[-1] = "  1.00"
                            this_occupancy[-1] = "  1.00"
                    except:
                        occupancy.append("  0.00")
                        this_occupancy = "  0.00"
                    try:
                        beta.append(string.strip(lin[60:66]))		#	61-66		temperature factor
                        this_beta = string.strip(lin[60:66])		#	61-66		temperature factor
                        if(beta[-1] == ''):
                            beta[-1] = "  0.00"
                    except:
                        beta.append("  0.00")
                        this_beta = "  0.00"
                    try:
                        segname.append(string.strip(lin[72:76]))	#	73-76		segment identifier
                        this_segname = string.strip(lin[72:76])	#	73-76		segment identifier
                        if(segname[-1] == '' and this_chain !=''):
                            segname[-1] = this_chain
                            this_segname = this_chain
                    except:
                        this_segname = ""
                        segname.append("")
                    try:
                        element.append(string.strip(lin[76:78]))	#	77-78		element symbol
                        if(element[-1] == ''):
                            element[-1] = "  "
                    except:
                        element.append("  ")
                    try:
                        charge.append(string.strip(lin[78:80]))		#	79-80		charge on the atom
                        if(charge[-1] == ''):
                            charge[-1] = "  "
                    except:
                        charge.append("  ")

                else:
                    occupancy.append(string.strip(lin[54:60]))      #	55-60		occupancy
                    this_occupancy =string.strip(lin[54:60])      #	55-60		occupancy
                    beta.append(string.strip(lin[60:66]))		#	61-66		temperature factor
                    this_beta = string.strip(lin[60:66])		#	61-66		temperature factor
                    segname.append(string.strip(lin[72:76]))	#	73-76		segment identifier
                    this_segname = string.strip(lin[72:76])	#	73-76		segment identifier
                    element.append(string.strip(lin[76:78]))	#	77-78		element symbol
                    charge.append(string.strip(lin[78:80]))		#	79-80		charge on the atom

                if(this_name not in unique_names): unique_names.append(this_name)
                if(this_resname not in unique_resnames): unique_resnames.append(this_resname)
                if(this_resid not in unique_resids): unique_resids.append(this_resid)
                if(this_chain not in unique_chains): unique_chains.append(this_chain)
                if(this_segname not in unique_segnames): unique_segnames.append(this_segname)
                if(this_occupancy not in unique_occupancies): unique_occupancies.append(this_occupancy)
                if(this_beta not in unique_betas): unique_betas.append(this_beta)
	
                this_resname=(string.strip(lin[17:21]))
                if this_resname in protein_resnames:
                    moltype.append('protein')
                    this_moltype = 'protein'
                elif this_resname in rna_resnames:	
                    moltype.append('rna')
                    this_moltype = 'rna'
                elif this_resname in dna_resnames:	
                    moltype.append('dna')
                    this_moltype = 'dna'
                elif this_resname in water_resnames:	
                    moltype.append('water')
                    this_moltype = 'water'
                else:
                    moltype.append('other')
                    this_moltype = 'other'
				
                if(this_moltype not in unique_moltypes): unique_moltypes.append(this_moltype)
					
                if(true_index == num_atoms):				
                    if(printme): print 'finished reading frame = ',this_frame
                    index=numpy.array(index,numpy.int)
                    original_index=numpy.array(original_index,numpy.int)
                    resid=numpy.array(resid,numpy.int)
                    original_resid=numpy.array(original_resid,numpy.int)

                    x=numpy.array(x,numpy.float32)
                    y=numpy.array(y,numpy.float32)
                    z=numpy.array(z,numpy.float32)
                    coor[0,:,0]=x.astype(numpy.float) ; coor[0,:,1]=y.astype(numpy.float) ; coor[0,:,2]=z.astype(numpy.float)
                    true_index = 0
                    x=[] ; y=[] ; z=[]
                    if(this_frame == 1):
                        self._atom=atom ; self._index=index  ; self._original_index = original_index ; self._name=name ; self._loc=loc ; self._resname=resname ; self._residue_flag = residue_flag
                        self._chain=chain ; self._resid=resid ; self._rescode=rescode ; self._original_resid=original_resid
                        self._occupancy=occupancy ; self._beta=beta ; self._segname=segname ; self._element=element
                        self._charge=charge ; self._moltype=moltype
                        
                        self._number_of_names = len(unique_names) ; self._names = unique_names
                        self._number_of_resnames = len(unique_resnames) ; self._resnames = unique_resnames
                        self._number_of_resids = len(unique_resids) ; self._resids = unique_resids
                        self._number_of_chains = len(unique_chains) ; self._chains = unique_chains
                        self._number_of_segnames = len(unique_segnames) ; self._segnames = unique_segnames
                        self._number_of_occupancies = len(unique_occupancies) ; self._occupancies = unique_occupancies
                        self._number_of_betas = len(unique_betas) ; self._betas = unique_betas
                        self._number_of_moltypes = len(unique_moltypes) ; self._moltypes = unique_moltypes
                        if(fastread):
                            self._natoms=len(index)
                            break

                    this_frame += 1
            elif((record_name != 'ATOM' or record_name != 'HETATM' and record_name != 'CONECT') and this_frame == 1):
                header.append(lin)
	
            elif((record_name == 'ATOM' or record_name == 'HETATM') and this_frame > 1):
                true_index += 1
                x.append(lin[30:38])				#	31-38		Real(8.3) X: angstroms	
                y.append(lin[38:46])				#	39-46		Real(8.3) Y: angstroms	
                z.append(lin[46:54])				#	47-54		Real(8.3) Z: angstroms	
		
                if(true_index == num_atoms):
                    if(printme): print 'finished reading frame = ',this_frame
                    x=numpy.array(x,numpy.float32)
                    y=numpy.array(y,numpy.float32)
                    z=numpy.array(z,numpy.float32)
                    coor[this_frame-1,:,0]=x ; coor[this_frame-1,:,1]=y ; coor[this_frame-1,:,2]=z
                    true_index = 0
                    this_frame += 1		
                    x=[] ; y=[] ; z=[]

            elif ((record_name == 'CONECT') and pdbscan):
                # Format of CONECT line:
                # Record name CONECT followed by list of atom indexes in 6 
                # character columns. First is the base atom, following atoms 
                # are those connected to it. 
                # Input line is filtered to ignore blank columns.
                ndxs = [int(lin[i:i+5]) for i in range(6, len(lin), 5) if lin[i:i+5].strip()]
                conect[ndxs[0]] = ndxs[1:]

        if(((this_frame - 1) != num_frames) and not fastread):
            if(printme): print '>>> WARNING: pdb file had ',num_frames,' sasio read ',this_frame-1,' frames'
            if(printme): print '>>> WARNING: pdb file had ',num_frames,' sasio read ',this_frame-1,' frames'
            if(printme): print '>>> WARNING: pdb file had ',num_frames,' sasio read ',this_frame-1,' frames'

        self._coor=numpy.array(coor)

        if 'check_zero_coor' in kwargs:
            self.check_for_all_zero_columns(self._coor)
			
        unique_elements = self.element_filter()

        self._number_of_elements = len(unique_elements) ; self._elements = unique_elements

        self._natoms=len(index)

        total_mass = self.calcmass()

        error = []

        if 'saspdbrx_topology' in kwargs:
            if kwargs['saspdbrx_topology']:	
                error = self.check_charmm_atomic_order_reorganize()
                return error

        self._header = header
        self._conect = conect

        return 

    def create_conect_pdb_lines(self):
        """
            Output stored conect information in PDB record format.
            Indices are converted from those read in - original_index - to 
            those used in output.
            Format of CONECT line:
            Record name CONECT followed by list of atom indexes in 6 character 
            columns. First is the base atom, following atoms are those 
            connected to it.
        """

        original_conect = self.conect()
        original_indexs = self.original_index()
        indexs = self.index()
            
        # Convert from original to current indexing
        convert_index = dict(zip(original_indexs, indexs))

        new_conect = {}

        for base, linked in original_conect.iteritems():

            new_base = convert_index[base]
            new_linked = []

            for index in linked:
                new_linked.append(convert_index[index])

                new_conect[new_base] = new_linked

        # Output should be in the order of the atom indexes
        ordered_keys = sorted(new_conect.keys())

        conect_lines = []

        for base in ordered_keys:
            ndxs = ''.join([str(i).rjust(5) for i in new_conect[base]])
            conect_lines.append('CONECT' + str(base).rjust(5) + ndxs)

        return conect_lines
