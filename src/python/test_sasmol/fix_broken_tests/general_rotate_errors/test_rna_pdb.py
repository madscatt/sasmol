import sasmol.sasmol as sasmol
import numpy

DataPath = '../../data/pdb_common/' 

mol = sasmol.SasMol(0)
mol.read_pdb(DataPath+'rna.pdb')
frame = 0
theta=numpy.pi/3.0

mol.general_axis_rotate(frame,theta,0.2,1.3,-3.5)
result_coor = mol.coor()
result_com  = mol.calccom(0)

print 'result_coor = ', result_coor
print 'result_com = ', result_com

