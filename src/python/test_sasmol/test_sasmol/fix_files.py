import glob,sys,os

files = glob.glob('*.py')

st1 = "SASSIE: Copyright (C) 2011 Joseph E. Curtis, Ph.D."
st2 = "Core-Testing: Copyright (C) 2011 Hailiang Zhang, Ph.D."
st3 = "from sassie.core_testing.util import env, util"
st4 = "from sassie.sasmol import sasmol, sasop, sascalc"
st5 = "import os; DataPath = os.path.dirname(os.path.realpath(__file__))+'/../../data/sasmol/sasmol/'"

s1c = s2c = s3c = s4c = s5c = 0
f = 0 ; l = 0
for my_file in files:
    if my_file == '__init__.py' or my_file == "fix_files.py":
        continue
    this_file = open(my_file,'r').readlines()
    print len(this_file) 
    for line in this_file:
        if f == 0 and l < 10:
            print line
        if line[:10] == st1[:10]:
               s1c += 1 
        elif line[:10] == st2[:10]:
               s2c += 1 
        elif line == st3[:10]:
               s3c += 1 
        elif line == st4:    
               s4c += 1 
        elif line == st5:    
               s5c += 1 
        l += 1
    sys.exit()

    f+=1
print 'number of files = ',len(files)
print 's1c = ',s1c
print 's2c = ',s2c
print 's3c = ',s3c
print 's4c = ',s4c
print 's5c = ',s5c





