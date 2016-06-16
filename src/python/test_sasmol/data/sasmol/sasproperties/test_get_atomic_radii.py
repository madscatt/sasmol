
import sasmol.sasmol as sasmol				
import numpy as np


#   setup the samol object
mol = sasmol.SasMol(0)
out_filename = "vdw.txt"
f = open(out_filename, 'w')

#   get the atomic radii

atomic=sasmol.sasproperties.Atomic()
#R_VDW = atomic.van_der_Waals_radii(keep_lower_case=True)
R_VDW = atomic.van_der_Waals_radii()


for key, item in R_VDW.items():	
    print key, item		
    f.write('%s' '%5s' '\n' %(key, item))

f.close()


