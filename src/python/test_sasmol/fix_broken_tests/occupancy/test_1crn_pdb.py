import sasmol.sasmol as sasmol
import numpy

DataPath = '../../data/pdb_common/' 

mol = sasmol.SasMol(0)
mol.read_pdb(DataPath+'1CRN.pdb')

occupancy = mol.occupancy()

print 'occupancy = ', occupancy

