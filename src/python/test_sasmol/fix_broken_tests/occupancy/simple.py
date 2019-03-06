import sasmol.sasmol as sasmol
import string
import numpy
import locale

DataPath = '../../data/pdb_common/' 

mol = sasmol.SasMol(0)
mol.read_pdb(DataPath+'1CRN.pdb')

filename = DataPath+'1CRN.pdb'

infile=open(filename,'r').readlines()

keep_going = True

while keep_going:
    
    for i in range(len(infile)):
        lin=infile[i]
        record_name = string.strip(lin[0:6])
        if((record_name[:4] == 'ATOM' or record_name == 'HETATM')):
            print 'lin = ', lin
            print 'lin[54:60] = ', lin[54:60] 
            string_occupancy = lin[54:60] 
            float_occupancy = locale.atof(lin[54:60]) 
            print 'float_occupancy = ', float_occupancy
            break

    keep_going = False


print 'mol.occupancy() = ', mol.occupancy()
print 'mol.occupancies() = ', mol.occupancies()
#print 'mol.beta() = ', mol.beta()
print 'mol.betas() = ', mol.betas()
#print 'mol.resid() = ', mol.resid()
#print 'mol.resids() = ', mol.resids()

if 11.81 in mol.betas():
    print 'yup it is there'
    print mol.betas().index(11.81)

if 13.79 in mol.betas():
    print 'yup it is there'
    print mol.betas().index(13.79)

resid = mol.resid()
unique_resids = mol.resids()
print 'unique_resids.index(resid[0]) = ', unique_resids.index(resid[i])

occupancy = mol.occupancy()
unique_occupancies = mol.occupancies()
print 'unique_occupancies.index(occupancy[0]) = ', unique_occupancies.index(occupancy[i])

beta = mol.beta()
nbeta = mol.beta().tolist()

unique_betas = mol.betas()
print 'beta[0] = ', beta[0]
print 'type(beta[0]) = ', type(beta[0])

print 'type(beta[0]) = ', type(beta[0])

#print [ "{:0.2f}".format(x) for x in a ]
nbeta = [ float("{:6.2f}".format(x)) for x in beta ]

print 'unique_betas.index(nbeta[0]) = ', unique_betas.index(nbeta[0])

#mol.initialize_children()
