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
import sasio
import sascalc
import sasop
import sassubset
import sasproperties
import saspdbrx
import sasview

#	SASMOL
#
#	12/4/2009	--	initial coding			:	jc
#	12/10/2009	--	doc strings 			:	jc
#	01/11/2010	--	new design pattern		:	jc
#
#	 1         2         3         4         5         6         7
#LC4567890123456789012345678901234567890123456789012345678901234567890123456789
#								       *      **
'''
	Sasmol is the main module that contains the base classes that 
	describe the objects in the system.  An instance of the system 
	can be composed of SasAtm and/or SasSol. 

	SasSys is the base class, it can accept atomic or solid types.
	
	For purely atomic based systems, a set of utilities are provided
	to hold information, to calculate properties, and to manipulate
	the structures in space.  In addition, basic file input / output
	is provided in the SasAtm class (and inherited up the classes
	as dictated by the calling code).

'''

class SasAtm(sasio.Files,sascalc.Prop,sasop.Move,sassubset.Mask,sasproperties.Atomic,saspdbrx.Topology,sasview.View):

    '''
	SasAtm is the base class to build and deal with atomistic systems.
	The various methods described herein are separated in regards to
	their function.  The class inherits file input/output and 
	manipulation abilities from the sasio, sascalc, and sasop modules.

    '''

    ###	OPEN	If you load an atom it doesn't check if it IS an atom 

    def __init__(self,id=None,filename=None):

###	OPEN	Test the init of ancestor init functions

        self._id=id
        self._totalmass=0.0	
        self._natoms=0
        self._mass=None
        self._coor=None
        self._com=None

    def setId(self,newValue):
        self._id = newValue

    def id(self):
        return self._id	

    def maketop(self):
        pass

    def setFilename(self,newValue):
        self._filename = newValue
	
    def filename(self):
        return self._filename
	
    def setType(self,newValue):
        self._type = newValue
	
    def type(self):
        return self._type

    def setCharmm_type(self,newValue):
        self._charmm_type = newValue
	
    def charmm_type(self):
        return self._charmm_type

	# file stuff

    def load(self):
        pass

    def delete(self):
        pass

    def save(self):
        pass

    def delframes(self):
        pass

	# properties

    def number_of_frames(self):
        return self._coor[:,0,0].shape[0]
	
    def setNumber_of_frames(self,newValue):
        self._number_of_frames = newValue
	
    def atom(self):
        return self._atom
	
    def setAtom(self,newValue):
        self._atom = newValue
	
    def index(self):
        return self._index
	
    def setIndex(self,newValue):
        self._index = newValue

    def header(self):
        return self._header
	
    def setHeader(self,newValue):
        self._header = newValue

    def original_index(self):
        return self._original_index
	
    def setOriginal_index(self,newValue):
        self._original_index = newValue

    def original_resid(self):
        return self._original_resid
	
    def setOriginal_resid(self,newValue):
        self._original_resid = newValue

    def conect(self):
        return self._conect
	
    def setConect(self,newValue):
        self._conect = newValue
	
    def residue_flag(self):
        try:
            return self._residue_flag
        except:
            self._residue_flag = False
            return self._residue_flag
	
    def setResidue_flag(self,newValue):
        self._residue_flag = newValue
	
    def name(self):
        return self._name
	
    def setName(self,newValue):
        self._name = newValue

    def loc(self):
        return self._loc
	
    def setLoc(self,newValue):
        self._loc = newValue
	
    def resname(self):
        return self._resname
	
    def setResname(self,newValue):
        self._resname = newValue
	
    def chain(self):
        return self._chain
	
    def setChain(self,newValue):
        self._chain = newValue
	
    def resid(self):
        return self._resid
	
    def setResid(self,newValue):
        self._resid = newValue

    def coor(self):
        return self._coor
	
    def setCoor(self,newValue):
        self._coor=newValue

    def rescode(self):
        return self._rescode
	
    def setRescode(self,newValue):
        self._rescode = newValue
	
    def occupancy(self):
        return self._occupancy
	
    def setOccupancy(self,newValue):
        self._occupancy = newValue
	
    def beta(self):
        return self._beta
	
    def setBeta(self,newValue):
        self._beta = newValue
	
    def segname(self):
        return self._segname
	
    def setSegname(self,newValue):
        self._segname = newValue
	
    def element(self):
        return self._element
	
    def setElement(self,newValue):
        self._element = newValue
	
    def charge(self):
        return self._charge	

    def setCharge(self,newValue):
        self._charge = newValue

    def atom_charge(self):
        return self._atom_charge	

    def setAtom_charge(self,newValue):
        self._atom_charge = newValue

    def atom_vdw(self):
        return self._atom_vdw	

    def setAtom_vdw(self,newValue):
        self._atom_vdw = newValue

    def residue_charge(self):
        return self._residue_charge	

    def setResidue_charge(self,newValue):
        self._residue_charge = newValue

    def energy(self):
        return self._energy

    def setEnergy(self,newValue):
        self._energy = newValue

    def formula(self):
        return self._formula

    def setFormula(self,newValue):
        self._formula = newValue

    def mass(self):
        return self._mass	

    def setMass(self,newValue):
        self._mass = newValue

    def totalmass(self):
        if(self._totalmass==None ):
            self._totalmass=sascalc.Prop.calcmass(self)
        return self._totalmass	

    def setTotalmass(self,newValue):
        self._totalmass = newValue

    def unitcell(self):
        return self._unitcell

    def setUnitcell(self,newValue):
        self._unitcell = newValue

    def com(self):
        return self._com
	
    def setCom(self,newValue):
        self._com = newValue

    def natoms(self):
        return self._natoms

    def setNatoms(self,newValue):
        self._natoms = newValue

    def rg(self):
        return self._rg

    def setRg(self,newValue):
        self._rg = newValue

    def pmi(self):
        return self._pmi

    def setPmi(self,newValue):
        self._pmi = newValue

    def minimum(self):
        return self._minimum

    def setMinimum(self,newValue):
        self._minimum = newValue
	
    def maximum(self):
        return self._maximum

    def setMaximum(self,newValue):
        self._maximum = newValue
	
    def shape(self):
        return self._shape

    def setShape(self,newValue):
        self._shape=newValue
	
    def moltype(self):
        return self._moltype

    def setMoltype(self,newValue):
        self._moltype=newValue

    def number_of_names(self):
        return self._number_of_names
	
    def setNumber_of_names(self,newValue):
        self._number_of_names = newValue

    def number_of_resnames(self):
        return self._number_of_resnames
	
    def setNumber_of_resnames(self,newValue):
        self._number_of_resnames = newValue

    def number_of_resids(self):
        return self._number_of_resids

    def setNumber_of_resids(self,newValue):
        self._number_of_resids = newValue

    def number_of_chains(self):
        return self._number_of_chains
	
    def setNumber_of_chains(self,newValue):
        self._number_of_chains = newValue

    def number_of_segnames(self):
        return self._number_of_segnames
	
    def setNumber_of_segnames(self,newValue):
        self._number_of_segnames = newValue

    def number_of_occupancies(self):
        return self._number_of_occupancies
	
    def setNumber_of_occupancies(self,newValue):
        self._number_of_occupancies = newValue

    def number_of_betas(self):
        return self._number_of_betas
	
    def setNumber_of_betas(self,newValue):
        self._number_of_betas = newValue

    def number_of_moltypes(self):
        return self._number_of_moltypes
	
    def setNumber_of_moltypes(self,newValue):
        self._number_of_moltypes = newValue

    def number_of_elements(self):
        return self._number_of_elements
	
    def setNumber_of_elements(self,newValue):
        self._number_of_elements = newValue

    def names(self):
        return self._names
	
    def setNames(self,newValue):
        self._names = newValue

    def resnames(self):
        return self._resnames
	
    def setResnames(self,newValue):
        self._resnames = newValue

    def resids(self):
        return self._resids
	
    def setResids(self,newValue):
        self._resids = newValue

    def chains(self):
        return self._chains
	
    def setChains(self,newValue):
        self._chains = newValue

    def segnames(self):
        return self._segnames
	
    def setSegnames(self,newValue):
        self._segnames = newValue

    def occupancies(self):
        return self._occupancies
	
    def setOccupancies(self,newValue):
        self._occupancies = newValue

    def betas(self):
        return self._betas
	
    def setBetas(self,newValue):
        self._betas = newValue

    def moltypes(self):
        return self._moltypes
	
    def setMoltypes(self,newValue):
        self._moltypes = newValue

    def elements(self):
        return self._elements
	
    def setElements(self,newValue):
        self._elements = newValue

    def names_mask(self):
        return self._names_mask
	
    def setNames_mask(self,newValue):
        self._names_mask = newValue

    def resnames_mask(self):
        return self._resnames_mask
	
    def setResnames_mask(self,newValue):
        self._resnames_mask = newValue

    def resids_mask(self):
        return self._resids_mask
	
    def setResids_mask(self,newValue):
        self._resids_mask = newValue

    def chains_mask(self):
        return self._chains_mask
	
    def setChains_mask(self,newValue):
        self._chains_mask = newValue

    def occupancies_mask(self):
        return self._occupancies_mask
	
    def setOccupancies_mask(self,newValue):
        self._occupancies_mask = newValue

    def betas_mask(self):
        return self._betas_mask
    
    def setBetas_mask(self,newValue):
        self._betas_mask = newValue

    def elements_mask(self):
        return self._elements_mask
	
    def setElements_mask(self,newValue):
        self._elements_mask = newValue

    def segnames_mask(self):
        return self._segnames_mask
	
    def setSegnames_mask(self,newValue):
        self._segnames_mask = newValue

    def effective_charge(self):
        return self._effective_charge
	
    def setEffective_charge(self,newValue):
        self._effective_charge = newValue

    def long_range_potential(self):
        return self._long_range_potential
	
    def setLong_range_potential(self,newValue):
        self._long_range_potential = newValue

    def short_range_potential(self):
        return self._short_range_potential
	
    def setShort_range_potential(self,newValue):
        self._short_range_potential = newValue

    def effective_van_der_waals(self):
        return self._effective_van_der_waals
	
    def setEffective_van_der_waals(self,newValue):
        self._effective_van_der_waals = newValue

    def fasta(self):
        return self._fasta
	
    def setFasta(self,newValue):
        self._fasta = newValue

    def one_letter_resname(self):
        return self._one_letter_resname
	
    def setOne_letter_resname(self,newValue):
        self._one_letter_resname = newValue

	### OPEN	Add an iterator to count over the atoms ????

class SasMol(SasAtm):

    '''
    SasMol is a class that is used to describe molecules. It inherits
    all of the attributes from SasAtm.  An example of a molecule is
    a single protein, a single nucleic acid strand.

    '''

###	OPEN	If you load a molecule it doesn't check to see if it IS an molecule

    def __init__(self,id):
        SasAtm.__init__(self,id)
        self._name='Mol_None'

    def molcharge(self):
        pass

    def corr(self):
        pass

class SasAss(SasMol):

    '''
    SasAss is a class that is used to describe atomic assemblies
    It inherits all of the attributes from SasMol and SasAtm.  An
    example of an assembly is two protein molecules, a box of water,
    a lipid membrane bilayer.  An assembly can be built by the
    user or it can be read in a single file.
    '''
###  	DONE    Ans: Yes. If SasAss can be loaded directly (i.e. with SasMol) does it inherit from SasAtm?

###	OPEN	If you load an assembly it doesn't check to see if it IS an assembly

    def __init__(self,id):
        SasMol.__init__(self,id)
        self._name='Ass_None'

class SasSol(object):

    '''
    SasSol is a class used to describe a geometric object in the
    system that is not easily desribed by an atomic representation.
    An example of a solid is a plane, a sphere, an irregular solid
    shape.
    '''

    def __init__(self,id):
        self._name='Sol_None'

class SasHybrid(SasMol):

    '''
    SasHybrid is a class used to describe a mixed atomic/geometric 
    system that is not easily desribed by an atomic representation
    alone.
    '''
    
    def __init__(self,id):
        self._name='Hybrid_None'

class SasSys():

    '''
    SasSys is a class that is used to describe a complete system.
    Deletion of an instance of SasSys removes all information and
    memory usage related to that object.

    Usage:
        import sasmol
        id=0
        a=sasmol.SasSys(id)   		# creates an atomic object by default (see comment on next line)
        a=sasmol.SasSys(id,'atomic') 	# creates an atomic object; system will only contain atomic objects
        a=sasmol.SasSys(id,'solid') 	# creates a solid type; system will ONLY contain solid objects
        a=sasmol.SasSys(id,'hybrid') 	# creates a mixed atomic and solid system

    '''

    def __init__(self,id,systype='atomic'):
        self._id=id
        if(systype=='atomic' or systype==None):
            self._systype='atomic'
        elif(systype=='solid'):
            self._systype='solid'
        else:
            pass
###	OPEN	Exception for systype error

        self._objectarray=[]

    def id(self):
        return self._id	

    def systype(self):
        return self._systype

    def add_object(self,addtype='atomic'):
        if(addtype=='atomic'):
            thisid=len(self._objectarray)
            object=SasAtm(thisid)
            self._objectarray.append([thisid,object])
            return object
        elif(addtype=='solid'):
            thisid=len(self._objectarray)
            object=SasSol(thisid)
            self._objectarray.append([thisid,object])
            return object
        elif(addtype=='hybrid'):
            thisid=len(self._objectarray)
            object=SasHybrid(thisid)
            self._objectarray.append([thisid,object])
            return object
#
###	OPEN	Add solid objects here.
#
        else:
            print '>>> error: need to specify addtype == atom, molecule, or assembly'
            print '>>> no objects created'
#
###	OPEN	Error check for add_object "atomic" input parameters
#
###	RUSH	Need to figure out how to manipulate by IDs, object names, and maybe iterators
#
###	RUSH	Can NOT remove objects using a.remove_object yet... "del object name" works
#
    def remove_object(self,object=None,id=None):
        print 'objectarray = ',self._objectarray
        del object
        for i in range(len(self._objectarray)):
            thisid=self._objectarray[i][0]	
            if thisid==id:
                self._objectarray[i][1]=None
                break
            #else:
            #	print 'Incorrect object id and/or object specified'

        print 'objectarray = ',self._objectarray

