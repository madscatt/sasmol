
import sasmol.sasmol as sasmol				


#   setup the samol object
mol = sasmol.SasMol(0)
out_filename = "dna_sl.txt"
f = open(out_filename, 'w')


atomic=sasmol.sasproperties.Atomic()

dsl = atomic.dna_sl()


for key, item in dsl.items():	
    print key, item		
    f.write('%3s' ' ' '%s' '\n' %(key, item))

f.close()


