#include <stdlib.h>
#include <stdio.h>
#include <string.h>
/*
    SASSIE  Copyright (C) 2011 Joseph E. Curtis
    This program comes with ABSOLUTELY NO WARRANTY; 
    This is free software, and you are welcome to redistribute it under certain
    conditions; see http://www.gnu.org/licenses/gpl-3.0.html for details.
*/
void get_mask_array(long long *aptr, int nflexible, int natoms, char **name,long long *resid,long long *flexible_residues,int nresidues, int mtype){

	int i, j, q0, fr, value, count ;
	long long value_array[natoms] ;

	long long (*farray)[natoms] = (long long(*)[natoms])aptr ;

	for(i=0;i<natoms;i++){
		value_array[i]=0 ;
	}

        for(i=0;i<nflexible;i++) {
		for(j=0;j<natoms;j++){
			farray[i][j]=0;
		}
	}

	count=0 ;

	if(mtype == 0) {
		for(fr=0 ; fr<nflexible ; fr++){
			q0 = flexible_residues[fr] ;
			for(i=0 ; i<natoms ; i++){
				value=0 ;
				if (resid[i]<q0-1 || resid[i]>q0+1)
                		{
                    			count++;
                    			value_array[i]=0;
                    			continue;
                		}
				if(resid[i] == q0-1 && (strcmp(name[i],"C")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"N")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"CA")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"C")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i],"N")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else{
					value_array[i]=value ;
				}
			}
			for(i=0 ; i<natoms ; i++){
				farray[fr][i]=value_array[i] ;		
			} 
		} 
	}
	else if(mtype == 1){
		for(fr=0 ; fr<nflexible ; fr++){
			q0 = flexible_residues[fr] ;
			for(i=0 ; i<natoms ; i++){
				value=0 ;
				if (resid[i]<q0-1 || resid[i]>q0+1)
                		{
                    			count++;
                    			value_array[i]=0;
                    			continue;
                		}
				if(resid[i] == q0-1 && (strcmp(name[i],"O3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"P")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"O5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"C5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"C4'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"C3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0 && (strcmp(name[i],"O3'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i],"P")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else if(resid[i] == q0+1 && (strcmp(name[i],"O5'")==0)){
					value=1 ; 
					count++;
					value_array[i]=value ;
				}
				else{
					value_array[i]=value ;
				}
			}
			for(i=0 ; i<natoms ; i++){
				farray[fr][i]=value_array[i] ;		
			} 
		} 
	}

	return ;

} ; /* end of get_mask_array */

/* for testing */

/*
void main() {

	int i,fr,nresidues,natoms,nflexible ;
	int resid[30], flexible_residues[2] ;
	char *name[30] ;
	int *aptr ;
	
	nflexible=2 ;	
	nresidues=6 ;
	natoms=30 ;

	int ** farray = (int **)malloc(nflexible*sizeof(int*));
	memset(farray,0,nflexible*sizeof(int*));
	for(i=0;i<nflexible;i++) {
		farray[i] = (int*)malloc(natoms*sizeof(int));
	}

	aptr= &farray[0][0] ;

	flexible_residues[0]=4 ;
	flexible_residues[1]=5 ;
	
	resid[0]=1 ; resid[1]=1 ; resid[2]=1 ; resid[3]=1 ; resid[4]=1 ;
	resid[5]=2 ; resid[6]=2 ; resid[7]=2 ; resid[8]=2 ; resid[9]=2 ;
	resid[10]=3 ; resid[11]=3 ; resid[12]=3 ; resid[13]=3 ; resid[14]=3 ;
	resid[15]=4 ; resid[16]=4 ; resid[17]=4 ; resid[18]=4 ; resid[19]=4 ;
	resid[20]=5 ; resid[21]=5 ; resid[22]=5 ; resid[23]=5 ; resid[24]=5 ;
	resid[25]=6 ; resid[26]=6 ; resid[27]=6 ; resid[28]=6 ; resid[29]=6 ;

	name[0]="N" ; name[1]="CA" ; name[2]="C" ; name[3]="HN" ; name[4]="HA" ;
	name[5]="F" ; name[6]="GA" ; name[7]="CA" ; name[8]="HN" ; name[9]="HA" ;
	name[10]="F" ; name[11]="CA" ; name[12]="NA" ; name[13]="HN" ; name[14]="HA" ;
	name[15]="N" ; name[16]="C" ; name[17]="N" ; name[18]="HN" ; name[19]="CA" ;
	name[20]="F" ; name[21]="CA" ; name[22]="NA" ; name[23]="HN" ; name[24]="HA" ;
	name[25]="N" ; name[26]="C" ; name[27]="N" ; name[28]="HN" ; name[29]="CA" ;
	
	get_mask_array(aptr,nflexible,natoms,name,resid,flexible_residues,nresidues);

	printf("back in main\n");
	for(fr=0;fr<nflexible;fr++){
		for(i=0;i<natoms;i++){
			if(farray[fr][i] != 0){
				printf("%i\t%i\t%i\n",fr,i,farray[fr][i]) ;
			}
		}
		printf("\n");
	}

	printf("\n");

} 
*/

/* end of main */
