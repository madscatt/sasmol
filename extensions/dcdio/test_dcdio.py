'''
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
'''
import sys,numpy,time
sys.path.append('./')
import dcdio
import sasmol

A=sasmol.SasMol(0)
A.read_pdb('min3.pdb')

natoms = A.natoms()

x=A.coor()[:,0]
y=A.coor()[:,1]
z=A.coor()[:,2]

#x=[1.0,2.0,3.0] ; y=[1.0,2.0,3.0] ; z=[1.0,2.0,3.0]

x=numpy.array(x,numpy.float32)
y=numpy.array(y,numpy.float32)
z=numpy.array(z,numpy.float32)

filename='c7.dcd'

fp=dcdio.open_dcd_write(filename)

nset=200; istart=1 ; nsavc=1 ; delta=1.0

headerresult=dcdio.write_dcdheader(fp,filename,natoms,nset,istart,nsavc,delta)

print 'writing '+str(nset)+' to disk'
start_time=time.time()

for blah in range(nset):
     print ".",
     sys.stdout.flush()

     x=x+5.0
     y=y+5.0
     z=z+5.0
     stepresult=dcdio.write_dcdstep(fp,x,y,z,blah)

end_time=time.time()

dt=end_time-start_time

print '\ntotal time = ',dt,' time per structure = ',dt/nset

dcdio.close_dcd_write(fp)


filename='200c.dcd'
filename='c7.dcd'

ifp=dcdio.open_dcd_read(filename)

nnatoms=0 ; nset=0 ; istart=0 ; nsavc=0 ; delta=0.0
namnf=0 ; freeindexes=[] ; reverseEndian=0 ; charmm=0

print 'nnatoms = ',nnatoms
print 'nset = ',nset
print 'freeindexes = ',freeindexes

readheaderresult,nnatoms,nset,istart,nsavc,delta,namnf,reverseEndian,charmm=dcdio.read_dcdheader(ifp)

print 'read header result = ',readheaderresult

print 'nnatoms = ',nnatoms
print 'nset = ',nset
print 'istart = ',istart
print 'nsavc = ',nsavc
print 'delta = ',delta
print 'namnf = ',namnf
print 'reverseEndian = ',reverseEndian
print 'charmm = ',charmm

x=numpy.zeros((nset,nnatoms),dtype=numpy.float32)
y=numpy.zeros((nset,nnatoms),dtype=numpy.float32)
z=numpy.zeros((nset,nnatoms),dtype=numpy.float32)

num_fixed=0 ; first=1
result=1

i=0

#try:
print 'reading dcd file'
start_time=time.time()
sum=0.0
for i in xrange(nset):
	print '.',
 	sys.stdout.flush()
 	read_start_time=time.time()
	tx=numpy.zeros(nnatoms,dtype=numpy.float32)
	ty=numpy.zeros(nnatoms,dtype=numpy.float32)
	tz=numpy.zeros(nnatoms,dtype=numpy.float32)

	result=dcdio.read_dcdstep(ifp,tx,ty,tz,num_fixed,i,reverseEndian,charmm)
 	read_end_time=time.time()
	
	sum+=read_end_time-read_start_time

	x[i][:]=tx ; y[i][:]=ty ; z[i][:]=tz
	
end_time=time.time()
dt=end_time-start_time

print '\nread total_time = ',sum,' time per structure = ',sum/nset
print 'total_time = ',dt,' time per structure = ',dt/nset
print 'ratio(total time) = ',dt/sum
print 'ratio(per structure) = ',(dt/nset)/(sum/nset)


#except:
#
#	print " I failed          :(     "

dcdio.close_dcd_read(ifp)

