        subroutine matrix_multiply(a,b,c,dim_a1,dim_a2,dim_b2)
        integer dim_a1,dim_a2,dim_b2
        double precision a(dim_a1,dim_a2)
        double precision b(dim_a2,dim_b2)
        double precision c(dim_a1,dim_b2)
        double precision cij

cf2py	intent(in) :: a,b,dim_a1,dim_a2,dim_b1,dim_b2
cf2py	intent(out):: c
cf2py	intent(hide):: i,j,k

        do 98,i=1,dim_a1
                do 99, j=1,dim_b2
                        c(i,j) = 0d0

  99         continue
  98    continue

        do 300,i=1,dim_a1
           do 200,j=1,dim_b2
                cij = 0d0
                do 100, k=1,dim_a2
                        cij = cij + a(i,k)*b(k,j)

  100           continue

                c(i,j) = cij
 
  200        continue

  300   continue

        end
