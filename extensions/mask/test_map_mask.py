'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import os,sys,numpy
import sasmol

def old_way():

	m1=sasmol.SasMol(0)

	m1.read_pdb('min3.pdb')

	moltype=m1.moltype()

	natoms=m1.natoms()

	nresidues=m1.resid()[-1] - m1.resid()[0]+1

	flexible_residues=numpy.arange(2,8).tolist()

	print 'getting rotation indices for molecule'

	residue_rotation_indices={}
	residue_rotation_mask={}
	
	print 'number of flexible residues = ',len(flexible_residues)
	
	for q0 in xrange(2,nresidues):
		if q0 in flexible_residues:
        		print q0,
			sys.stdout.flush()
			previous_amino_acid=q0-1
			this_amino_acid=q0
			next_amino_acid=q0+1
	
			refatm = '(resid[i] == '+str(previous_amino_acid)+' and name[i] == "C")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "N")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "CA")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "C")'
			refatm += ' or (resid[i] == '+str(next_amino_acid)+' and name[i] == "N")'
	
			error,mask = m1.get_subset_mask(refatm)  ; indices = m1.get_indices_from_mask(mask)
	
			residue_rotation_indices[q0] = indices.tolist()
			residue_rotation_mask[q0] = mask.tolist()

def list_way():

	m1=sasmol.SasMol(0)

	m1.read_pdb('min3.pdb')

	moltype=m1.moltype()

	natoms=m1.natoms()

	resid=m1.resid()
	name=m1.name()

	nresidues=m1.resid()[-1] - m1.resid()[0]+1

	flexible_residues=numpy.arange(2,8).tolist()

	print 'getting rotation indices for molecule'

	residue_rotation_indices={}
	residue_rotation_mask={}
	
	print 'number of flexible residues = ',len(flexible_residues)

	atomlist=numpy.arange(0,natoms).tolist()
	
	for q0 in xrange(2,nresidues):
		if q0 in flexible_residues:
        		print q0,
			sys.stdout.flush()
			previous_amino_acid=q0-1
			this_amino_acid=q0
			next_amino_acid=q0+1
	
			refatm = '(resid[i] == '+str(previous_amino_acid)+' and name[i] == "C")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "N")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "CA")'
			refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "C")'
			refatm += ' or (resid[i] == '+str(next_amino_acid)+' and name[i] == "N")'

			torf = [eval(refatm) for i in atomlist]

#			torf = [eval(refatm) for i in atomlist]
#			torf = map(eval(refatm), atomlist)
	
#			map(eval,refatm)

#			error,mask = m1.get_subset_mask(refatm)  ; indices = m1.get_indices_from_mask(mask)
	
#			residue_rotation_indices[q0] = indices.tolist()
#			residue_rotation_mask[q0] = mask.tolist()


class StringFunction_0:
#	def __init__(self,expression):
	def __init__(self,expression,**kwargs):
		self._f_compiled = compile(expression, '<string>', 'eval')
		self._var = kwargs.get('name','resid')
		self._prms = kwargs
	
	def set_parameters(self, **kwargs):
		self.prms.update(kwargs)

	def __call__(self,*args):
		for n, value in zip(self._var,args):
			self._prms[n] = value	
		return eval(self._f_compiled, globals(), self._prms)

def new_way():

	m1=sasmol.SasMol(0)

	m1.read_pdb('min3.pdb')

	moltype=m1.moltype()

	natoms=m1.natoms()

	aresid=m1.resid()
	aname=m1.name()

	nresidues=m1.resid()[-1] - m1.resid()[0]+1

	flexible_residues=numpy.arange(2,8).tolist()

	print 'getting rotation indices for molecule'

	residue_rotation_indices={}
	residue_rotation_mask={}
	
	print 'number of flexible residues = ',len(flexible_residues)

	atomlist=numpy.arange(0,natoms).tolist()

	for q0 in xrange(2,nresidues):
		if q0 in flexible_residues:
        		print q0,
			sys.stdout.flush()
			previous_amino_acid=q0-1
			this_amino_acid=q0
			next_amino_acid=q0+1
	
			#refatm = '(resid[i] == '+str(previous_amino_acid)+' and name[i] == "C")'
			#refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "N")'
			#refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "CA")'
			#refatm += ' or (resid[i] == '+str(this_amino_acid)+' and name[i] == "C")'
			#refatm += ' or (resid[i] == '+str(next_amino_acid)+' and name[i] == "N")'

			refatm = '(resid == '+str(previous_amino_acid)+' and name == "C")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "N")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "CA")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "C")'
			refatm += ' or (resid == '+str(next_amino_acid)+' and name == "N")'

			preliminary_mask_array=[]
			for i in xrange(natoms):
				name=aname[i] ; resid=aresid[i]
				ans=StringFunction_0(refatm,name=name,resid=resid)
				if(ans(name,resid)==1):
				#	print q0,i
					preliminary_mask_array.append(1)
				else:
					preliminary_mask_array.append(0)

#			torf=[StringFunction('refatm',independent_variable='i') for i in atomlist]
	
			mask_array = numpy.array(preliminary_mask_array, numpy.int32)

#			error,mask = m1.get_subset_mask(refatm)  ; indices = m1.get_indices_from_mask(mask)
	
#			residue_rotation_indices[q0] = indices.tolist()
			residue_rotation_mask[q0] = mask_array.tolist()

#	print 'residue_rotation_mask[2] = ',residue_rotation_mask[2]


class LamFunction:
	def _build_lambda(self):
		s = 'lambda ' + ', '.join(self._var)
		if self._prms:
			s += ', ' + ', '.join(['%s=%s' % (k, self._prms[k]) \
				for k in self._prms])
		s += ': ' + self._f
		self.__call__ = eval(s, self._globals)
		

def lam_way():
	
#	import StringFunction

	m1=sasmol.SasMol(0) ; m1.read_pdb('min3.pdb') ; natoms=m1.natoms()
	aresid=m1.resid() ; aname=m1.name()
	nresidues=m1.resid()[-1] - m1.resid()[0]+1
	flexible_residues=numpy.arange(2,8).tolist()
	print 'getting rotation indices for molecule'
	residue_rotation_indices={} ; residue_rotation_mask={}
	print 'number of flexible residues = ',len(flexible_residues)
	atomlist=numpy.arange(0,natoms).tolist()

	for q0 in xrange(2,nresidues):
		if q0 in flexible_residues:
        		print q0,
			sys.stdout.flush()
			previous_amino_acid=q0-1
			this_amino_acid=q0
			next_amino_acid=q0+1
	
			refatm = '(resid == '+str(previous_amino_acid)+' and name == "C")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "N")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "CA")'
			refatm += ' or (resid == '+str(this_amino_acid)+' and name == "C")'
			refatm += ' or (resid == '+str(next_amino_acid)+' and name == "N")'

			preliminary_mask_array=[]
			for i in xrange(natoms):
				name=aname[i] ; resid=aresid[i]
#				ans=StringFunction.StringFunction(refatm, independent_variables=('name','resid'))
				if(ans(name,resid)==1):
				#	print q0,i
					preliminary_mask_array.append(1)
				else:
					preliminary_mask_array.append(0)

			mask_array = numpy.array(preliminary_mask_array, numpy.int32)
			residue_rotation_mask[q0] = mask_array.tolist()

#	print 'residue_rotation_mask[2] = ',residue_rotation_mask[2]

def cee_way():

	m1=sasmol.SasMol(0) ; m1.read_pdb('min3.pdb') ; natoms=m1.natoms()
	resid=m1.resid() ; name=m1.name()
	nresidues=m1.resid()[-1] - m1.resid()[0]+1
	flexible_residues=numpy.arange(2,8).tolist()
#	print 'getting rotation indices for molecule'
	residue_rotation_indices={} ; residue_rotation_mask={}
#	print 'number of flexible residues = ',len(flexible_residues)
	nflexible=len(flexible_residues)
	atomlist=numpy.arange(0,natoms).tolist()
	
	#import cee_mask
	import mask
	
	farray=numpy.zeros((nflexible,natoms),numpy.int32)

	nresidues=int(nresidues) 

	mask.get_mask_array(farray,name,resid,flexible_residues,nresidues)

	for i in xrange(nflexible):
		for j in xrange(natoms):
			if(farray[i][j]!=0):
				print farray[i][j],

#void get_mask_array(int **farray, char **name,int *resid,int *flexible_residues,int nresidues,int natoms,int nflexible){

if __name__=='__main__': 

	from timeit import Timer
# 	print 'LAM WAY'
# 	t=Timer('lam_way()', 'from __main__ import lam_way')
# 	print '\n',t.timeit(number=1)
# 	print '\n\nNEW WAY'
# 	nt=Timer('new_way()', 'from __main__ import new_way')
# 	print '\n',nt.timeit(number=1)
 	print '\n\nOLD WAY'
 	t=Timer('old_way()', 'from __main__ import old_way')
 	old_time=t.timeit(number=9)
 	print '\n',old_time
	print '\n\nCEE WAY'
	t=Timer('cee_way()', 'from __main__ import cee_way')
	cee_time=t.timeit(number=9)
	print '\n',cee_time
 	print '\nratio = ',old_time/cee_time	
		
