#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readInput(char * fileName, double * & val, int * &rowind,int * & colptr, int & m, int&  n, int & nzl){
	int i,j;
	FILE * fp= fopen(fileName, "r");
	if(fp==NULL)
	{
		printf("file read error!\n");
	}
	/*first line */
	fscanf(fp,"%d",&m);
	fscanf(fp,"%d",&n);
	fscanf(fp,"%d",&nzl);
	val = new double[nzl];
	rowind = new int[m+1];
	colptr = new int[nzl];
	
	/* second line */
	for(i=0;i<m+1;i++){
		fscanf(fp,"%d",&rowind[i]);
	}
	
	/*third line */
	for(i=0;i<nzl;i++){
		fscanf(fp,"%d",&colptr[i]);
	}
	
	/* fourth line */
	for(i=0;i<nzl;i++){
		fscanf(fp,"%lf",&val[i]);
		if (val[i]<1.e-10 && val[i]>-1.e-10) 
			val[i]=0.0;
	}
	fclose(fp);
}
