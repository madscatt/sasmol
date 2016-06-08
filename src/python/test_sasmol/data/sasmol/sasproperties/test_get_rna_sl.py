
import sasmol.sasmol as sasmol				


#   setup the samol object
mol = sasmol.SasMol(0)
out_filename = "rna_sl.txt"
f = open(out_filename, 'w')


atomic=sasmol.sasproperties.Atomic()

rsl = atomic.rna_sl()


for key, item in rsl.items():	
    print key, item		
    f.write('%3s' ' ' '%s' '\n' %(key, item))

f.close()


