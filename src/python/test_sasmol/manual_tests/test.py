import os
import sasmol.sasmol as sasmol

pdbDataPath = os.path.join(os.path.dirname(os.path.realpath(__file__)),'..','data','pdb_common')+os.path.sep

m = sasmol.SasMol(0)

try:
    m.read_pdb("hiv1_gag.pdb")
except:
    m.read_pdb(pdbDataPath+"hiv1_gag.pdb")


