#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void readInput(char * fileName,double * &val, int * &rowind,int * &colptr,int & m, int&  n, int & nzl);
#ifdef IN_TEST
int main(int argc, char* argv[])
{
	double* val;
	int *rowind;
	int *colptr;
	unsigned int *peIndex;
	int m, n, nzl;
	char *fileName;
	if (argc>1) {
		fileName=argv[1];
	}
	else {
		printf("give file name\n");
		return 0;
	}


	readInput(fileName,val, rowind, colptr, peIndex, m,  n, nzl);
	return 0;
}
#endif
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
		if (val[i]<1.e-10 && val[i]>-1.e-10) {
			//			printf("val changed %f\n", val[i]);
			val[i]=0.0;
		}
	}
	fclose(fp);
#ifdef IN_TEST
	for(int j=0; j<n; j++)
	{
		printf("column %d: ",j);
		for (int i=colptr[j]; i<colptr[j+1]; i++) {
			printf("%d: %f ",rowind[i],val[i]);
			if (rowind[i]<j) {
				printf("Not Lower tirangular!\n");
				return;
			}
		}
		printf("\n");
	}
#endif
/*	bool changed=false;
	for (int j=0; j<n; j++) {
		do {
			changed=false;
			for (int i=colptr[j]; i<colptr[j+1]-1; i++) {
				if (rowind[i]>rowind[i+1]) {
					int tmp=rowind[i];
					double tmpV=val[i];
					rowind[i]=rowind[i+1];
					val[i] = val[i+1];
					rowind[i+1]=tmp;
					val[i+1]=tmpV;
					changed=true;
				}
			}
		} while(changed);
	}
	// last element with nonzero in that row
	int *last_ind_in_row = new int[m];
	for (int i=0; i<m; i++) {
		last_ind_in_row[i] = -1;
	}
	for(int j=0; j<n; j++)
	{
		for (int i=colptr[j]; i<colptr[j+1]; i++) {
			if (last_ind_in_row[rowind[i]]!=-1) {
				// next is j
				peIndex[last_ind_in_row[rowind[i]]] |= j;
				// this is dependent
				peIndex[i] = 0x80000000;
			}
			last_ind_in_row[rowind[i]] = i;
		}
	}
	delete[] last_ind_in_row; 
#ifdef IN_TEST
	for (int i=0; i<nzl; i++) {
		printf("%d ", peIndex[i]);
	}
	printf("\n");
	for (int i=0; i<nzl; i++) {
		printf("%x ", peIndex[i]);
	}
	printf("\n");
#endif
 */
}

void readInput(char * fileName, int * &rowind,int * & colptr, int & m, int&  n, int & nzl){
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
//	val = new double[nzl];
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
//	for(i=0;i<nzl;i++){
//		fscanf(fp,"%lf",&val[i]);
//		if (val[i]<1.e-10 && val[i]>-1.e-10) {
			//			printf("val changed %f\n", val[i]);
//			val[i]=0.0;
//		}
//	}
	fclose(fp);
}


