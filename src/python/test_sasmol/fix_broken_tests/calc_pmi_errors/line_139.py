import sasmol.sasmol as sasmol
import numpy

DataPath = '../../data/pdb_common/' 

mol = sasmol.SasMol(0)

mol.setCoor(numpy.array([[[1.0, 2.0, 3.0],[4.0, 5.0, 6.0],[7.0, 8.0, 9.0],[1.0, 3.0, 5.0],[2.0, 4.0, 6.0],[0.0, 2.0, 3.0]]],numpy.float64))
mol.setElement(['C','N','O','C','N','O'])


mol.setNatoms(len(mol.element()))
result = mol.calcpmi(0)
print 'result = ', result
result_eigenvalues = result[0]
result_eigenvectors = result[1].T
result_I = result[2]

print 'result_eigenvalues = ', result_eigenvalues
print 'result_eigenvectors = ', result_eigenvectors
print 'result_I = ', result_I

