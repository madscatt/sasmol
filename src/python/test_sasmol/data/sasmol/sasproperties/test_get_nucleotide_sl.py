
import sasmol.sasmol as sasmol				


#   setup the samol object
mol = sasmol.SasMol(0)
out_filename = "nucleotide_sl.txt"
f = open(out_filename, 'w')


atomic=sasmol.sasproperties.Atomic()

nsl = atomic.nucleotide_sl()


for key, item in nsl.items():	
    print key, item		
    f.write('%3s' ' ' '%s' '\n' %(key, item))

f.close()


