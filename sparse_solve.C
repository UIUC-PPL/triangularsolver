/*
 author: Ehsan Totoni, totoni2@illinois.edu
 December, 2011
 */

#include <stdio.h>
#include <vector>
using namespace std;
#include "sparse_solve.decl.h"

#define MIN_ENTRIES_PER_X 20
/*readonly*/ CProxy_Main mainProxy;
/*readonly*/ int nElements;
/*readonly*/ char fileName[250];

struct row_attr {
	int chare; // no. of chare waiting
	int row; // local row index in that chare
};

#include "MessagePool.h"
#include "ColumnsSolve.h"

// for each block chare keeps deps
struct chare_deps_str {
	int chare_no; // no. of chare
	row_attr* nextRow; // next of each row
	int size;
};

/*mainchare*/
class Main : public CBase_Main
{
	CProxy_ColumnsSolve arr;
	double start_time;
	int total_columns;
	int total_chares; // all chares after creation
	bool check;
	double* xVal;
public:
	Main(CkArgMsg* m)
	{
		check=false;
		//Process command-line arguments
		if(m->argc < 3) {
			CkPrintf("arguments: num_chares filename\n");
			CkExit();
		}
		nElements=atoi(m->argv[1]);
		strcpy(fileName, m->argv[2]);
		
		delete m;
		CkPrintf("reading file: %s\n",fileName);	
		FILE * fp= fopen(fileName, "r");
		if(fp==NULL)
		{
			CkPrintf("file read error!\n");
			CkExit();
		}
		int ncols, nzl;
		fscanf(fp,"%d",&ncols);
		fscanf(fp,"%d",&ncols);
		fscanf(fp,"%d",&nzl);		
		fclose(fp);

		total_columns = ncols;
		//Start the computation
		CkPrintf("Running ColumnsSolve on %d processors for %d elements\n", CkNumPes(),nElements);
		mainProxy = thisProxy;
		//	CProxy_BlockMap myMap=CProxy_BlockMap::ckNew();
		CProxy_RRMap myMap=CProxy_RRMap::ckNew();
		//Make a new array using that map
		CkArrayOptions opts(nElements);
		opts.setMap(myMap);
		arr = CProxy_ColumnsSolve::ckNew(opts);
		CkCallback *cb = new CkCallback(CkIndex_Main::reportIn(), mainProxy);
		arr.ckSetReductionClient(cb);
		setup_input(fileName);

	};
	void setup_input(char* fileName)
	{
		 double* data; // data in columns sparse compressed
		int* rowInd;  // row index of each element
		int *colsInd; // index of each column in data array, one more than columns for convenience
		// read input
		int m,n, nzl;
		readInput(fileName, data, rowInd, colsInd ,m, n, nzl);
		// number of local columns
		int num_loc_cols = (total_columns/nElements);
		int *lastRowInd = new int[m];
		memset(lastRowInd, 0, m*sizeof(int));
		
		double *tmpData = new double[nzl]; // too much memory!!!!
		int *tmpCol = new int[nzl];
		int *tmpRow = new int[m+1];
		bool *tmpDep = new bool[m];
		int chareNo = nElements;
		int *map_rows = new int[m];

		vector<chare_deps_str> chare_deps; // nexts for all block chares
		row_attr* prev_in_row = new row_attr[m]; // previous chare in that row
		for (int i=0; i<m; i++) {
			prev_in_row[i].chare = -1;
		}
		int offdiags = 0;
		int diagdeps = 0;
		double total_analysis_time =0;

		for (int i=0; i<nElements; i++) {
			int lastNoDiagChare = chareNo;
			// first column
			int startCol = i*num_loc_cols;
			int endCol = (i+1)*num_loc_cols;
			if (i==nElements-1) {
				endCol = total_columns;
			}
			// min entries: constant times number of x values
			int min_entries = MIN_ENTRIES_PER_X*(endCol-startCol);
			int indep_row_no=0;
			double start_analysis = CmiWallTimer();
			reorder(map_rows, indep_row_no, startCol, endCol, prev_in_row,
					rowInd, colsInd, data, tmpData, tmpRow, tmpCol, tmpDep, chare_deps, i, lastRowInd);
			total_analysis_time += CmiWallTimer()-start_analysis;

			// last row of next diagonal chare
			int nextChareLastRow = endCol + (endCol-startCol);
			if (i>=nElements-2) {
				nextChareLastRow = total_columns;
			}
			int curr_row = endCol;
			// get size of below diagonal corresponding to next diagonal Chare
			int belowSize = getBelowSize(rowInd, colsInd, endCol, endCol, nextChareLastRow, lastRowInd);
			int entries = tmpRow[endCol-startCol];
			int firstBelowChunkRows = 0; // number of rows of first chunk below diagonal
			int firstBelowChunkMaxCol = 0; // last needed column (X value) for this part
			int restBelowChunkRows = 0; // number of rows after first chunk below diagonal
			int restBelowChunkMaxCol = 0; // last needed column (X value) for this part
			if (belowSize<min_entries && belowSize!=0) {				
				addBelowRows(firstBelowChunkRows, firstBelowChunkMaxCol, entries, rowInd, colsInd, data, tmpRow,
							 tmpCol, tmpData, tmpDep, startCol, endCol, endCol, nextChareLastRow, 
							 lastRowInd, map_rows, chare_deps, prev_in_row, i);
				offdiags += belowSize;
				curr_row = nextChareLastRow;
				// rest of rows
				belowSize = getBelowSize(rowInd, colsInd, endCol, nextChareLastRow, total_columns, lastRowInd);
				if (belowSize<min_entries && belowSize!=0) {
					restBelowChunkRows = firstBelowChunkRows; // for rowInd calcs
					addBelowRows(restBelowChunkRows, restBelowChunkMaxCol, entries, rowInd, colsInd, data, tmpRow,
								 tmpCol, tmpData, tmpDep, startCol, endCol, nextChareLastRow, total_columns,
								 lastRowInd, map_rows, chare_deps, prev_in_row, i);
					restBelowChunkRows -= firstBelowChunkRows;
					offdiags += belowSize;
					curr_row = total_columns;
				}
				if (belowSize==0) {
					curr_row = total_columns;
				}
			}
			
			chare_deps_str new_chare; new_chare.size=firstBelowChunkRows+restBelowChunkRows; new_chare.chare_no=i;
			new_chare.nextRow=new row_attr[firstBelowChunkRows+restBelowChunkRows];
			chare_deps.push_back(new_chare);
			
			int my_rows = endCol-startCol+firstBelowChunkRows+restBelowChunkRows;
			InputMsg* msg = new ((my_rows+1), entries, entries, my_rows) InputMsg;
			memcpy(msg->data, tmpData, entries*sizeof(double));
			memcpy(msg->colInd, tmpCol, entries*sizeof(int));
			memcpy(msg->rowInd, tmpRow,(my_rows+1)*sizeof(int));
			memcpy(msg->dep, tmpDep,(my_rows)*sizeof(bool));
			msg->num_cols = endCol-startCol;
			msg->num_rows = my_rows;
			msg->start_col = startCol;
			msg->diag = true;
			msg->indep_row_no = indep_row_no;
			msg->first_below_rows = firstBelowChunkRows;
			msg->first_below_max_col = firstBelowChunkMaxCol;
			msg->rest_below_rows = restBelowChunkRows;
			msg->rest_below_max_col = restBelowChunkMaxCol;
			arr[i].get_input(msg);
			int dep_cols = endCol-startCol-indep_row_no;
		//	CkPrintf("nzls:%d deps:%d total:%d frac:%f first_max:%d rest_max:%d first_rows:%d rest_rows:%d\n", entries, dep_cols, 
		//			endCol-startCol, dep_cols/(double)(endCol-startCol), firstBelowChunkMaxCol, restBelowChunkMaxCol, firstBelowChunkRows
		//			 ,restBelowChunkRows);
			diagdeps += dep_cols;
			
			int tmp_curr_row =0;

			while (curr_row!=m) {
				int tmp_curr_row =0;
				entries=0;
				while (entries<min_entries) {
					
					tmpRow[tmp_curr_row] = entries;
					int j = rowInd[curr_row]+lastRowInd[curr_row];
					while (colsInd[j]<endCol) {
						tmpData[entries] = data[j];
						tmpCol[entries] = map_rows[colsInd[j]-startCol];
						j++;
						entries++;
					}
					// if empty go to next row
					if (entries == tmpRow[tmp_curr_row]) {
						curr_row++;
						if (curr_row==m) {
							break;
						}
						continue;
					}
					lastRowInd[curr_row] = j-rowInd[curr_row];
					tmpDep[tmp_curr_row] = false;
					// if row is not empty, this is the last
					// set dependency of previous block chare
					int last_chare = prev_in_row[curr_row].chare;
					if (last_chare!=-1) {
						for (int k=0; k<chare_deps.size(); k++) {
							if (chare_deps[k].chare_no==last_chare) {
								row_attr tmprattr; tmprattr.chare=chareNo; tmprattr.row=tmp_curr_row;
								chare_deps[k].nextRow[prev_in_row[curr_row].row] = tmprattr;
							}
						}
						tmpDep[tmp_curr_row] = true;
					}
					row_attr rtr; rtr.chare= chareNo; rtr.row=tmp_curr_row;
					prev_in_row[curr_row] = rtr;
					
					curr_row++;
					tmp_curr_row++;
					
					if (curr_row==m) {
						break;
					}
				}
				// empty lower diagonal
				if (entries==0) {
					break;
				}
				offdiags += entries;
				// CkPrintf("chare:%d nonzeros:%d\n",chareNo, entries);
				tmpRow[tmp_curr_row] = entries;
				chare_deps_str new_chare; new_chare.size=tmp_curr_row; new_chare.chare_no=chareNo;
				new_chare.nextRow=new row_attr[tmp_curr_row];
				chare_deps.push_back(new_chare);
				InputMsg* msg = new ((tmp_curr_row+1), entries, entries, tmp_curr_row) InputMsg;
				memcpy(msg->dep, tmpDep, tmp_curr_row*sizeof(bool));
				memcpy(msg->data, tmpData, entries*sizeof(double));
				memcpy(msg->colInd, tmpCol, entries*sizeof(int));
				memcpy(msg->rowInd, tmpRow,(tmp_curr_row+1)*sizeof(int));
				msg->num_cols = num_loc_cols;
				msg->num_rows = tmp_curr_row;
				msg->diag = false;
				arr[chareNo].insert();
				arr[chareNo++].get_input(msg);
			}
			// send section proxy to diag element
			// printf("max_res_col:%d total:%d\n", restBelowChunkMaxCol, endCol-startCol);
			CProxySection_ColumnsSolve nondiags = CProxySection_ColumnsSolve::ckNew(arr, lastNoDiagChare, chareNo-1, 1);
			bool empty_nondiags = false;
			if(lastNoDiagChare== chareNo)
				empty_nondiags = true;
			arr[i].get_section(nondiags, empty_nondiags);
		}
		for (int i=0; i<chare_deps.size(); i++) {
			DepsMsg* msg = new (chare_deps[i].size) DepsMsg;
			memcpy(msg->deps, chare_deps[i].nextRow, chare_deps[i].size*sizeof(row_attr));	
			arr[chare_deps[i].chare_no].get_deps(msg);
		}
		CkPrintf("analysis time:%f\n",total_analysis_time/nElements);
		CkPrintf("offdiags:%d out of %d nonzeros, fraction:%f\n",offdiags, nzl, offdiags/(double)nzl);
		CkPrintf("diagdeps:%d out of %d rows, fraction:%f\n",diagdeps, m, diagdeps/(double)m);
		total_chares = chareNo;
		arr.init();
		delete[] map_rows;
		delete[] lastRowInd;
		delete[] tmpRow;
		delete[] tmpData;
		delete[] tmpCol;
		delete[] tmpDep;
		delete[] data;
		delete[] rowInd;
		delete[] colsInd;
		delete[] prev_in_row;
	}

	void reportIn()
	{
		CkPrintf("All done in %f\n",CmiWallTimer()-start_time);
		CkExit();
	}
	void initDone(void)
	{
		start_time = CmiWallTimer();
		arr.start();
	};

	void reorder(int *map_rows, int& no_indeps, int startCol, int endCol, row_attr* prev_in_row, int* rowInd,
				int* colInd, double* data,double* tmpData,int* tmpRow,int* tmpCol,
				 bool* tmpDep, vector<chare_deps_str> &chare_deps, int chare_no, int *lastRowInd)
	{
		int size=endCol-startCol;
		// depend on messages
		bool *mark_dep = new bool[size];
		memset(mark_dep, 0, size*sizeof(bool));
		// done reordering
		bool *mark_done = new bool[size];
		memset(mark_done, 0, size*sizeof(bool));
		
		// find dependents
		for (int i=startCol; i<endCol; i++) {
			if (prev_in_row[i].chare!=-1) { // receives message
				mark_dep[i-startCol] = true;
				continue;
			}
			for (int j=rowInd[i+1]-2; j>=rowInd[i]+lastRowInd[i]; j--) {
				// depends on dependant row
				if (mark_dep[colInd[j]-startCol]) {
					mark_dep[i-startCol] = true;
					break;
				}
			}
		}
		// place non-deps
		int target_nzls = 0;
		int target_row = 0;
		tmpRow[0]=0;
		for (int i=endCol-1; i>=startCol; i--) {
			if (!mark_dep[i-startCol] && !mark_done[i-startCol]) {
				place(i, target_row, target_nzls, map_rows, startCol, rowInd,
					  colInd, data, tmpData, tmpRow, tmpCol, mark_done, lastRowInd);
			}
		}
		// number of independent rows
		no_indeps = target_row;
		
		// place dependents
		for (int i=startCol; i<endCol; i++) {
			if (mark_dep[i-startCol]) {
				map_rows[i-startCol] = target_row;
				for (int j=rowInd[i]+lastRowInd[i]; j<rowInd[i+1]; j++) {
					tmpData[target_nzls] = data[j];
					tmpCol[target_nzls] = map_rows[colInd[j]-startCol];
					target_nzls++;
				}
				tmpRow[target_row+1] = target_nzls;
				if (prev_in_row[i].chare!=-1) {
					tmpDep[target_row-no_indeps] = true;
					for (int k=0; k<chare_deps.size(); k++) {
						if (chare_deps[k].chare_no==prev_in_row[i].chare) {
							row_attr tmprattr; tmprattr.chare=chare_no; tmprattr.row=target_row;
							chare_deps[k].nextRow[prev_in_row[i].row] = tmprattr;
						}
					}
				} else {
					tmpDep[target_row-no_indeps] = false;
				}

				target_row++;
			}
		}
		
	}
	void place(int row, int& target_row, int& target_nzls, int *map_rows, int startCol, int* rowInd,
			   int* colInd, double* data,double *tmpData,int* tmpRow,int* tmpCol, bool* mark_done, int *lastRowInd)
	{
		// place dependencies
		int j=rowInd[row+1]-2;
		for (; j>=rowInd[row]+lastRowInd[row]; j--) {
			if (colInd[j]<startCol || colInd[j]==row) {
				continue;
			}
			assert(rowInd[row+1]-2>=rowInd[row]);
			assert(lastRowInd[row]>=0);
			assert(colInd[j]<row && colInd[j]>=startCol);
			if (!mark_done[colInd[j]-startCol]) {
				place(colInd[j], target_row, target_nzls, map_rows, startCol, rowInd,
					  colInd, data, tmpData, tmpRow, tmpCol, mark_done, lastRowInd);
			}
		}
		map_rows[row-startCol] = target_row;
		mark_done[row-startCol] = true;
		
		for (int j= rowInd[row]+lastRowInd[row]; j<rowInd[row+1]; j++) {
			tmpData[target_nzls] = data[j];
			tmpCol[target_nzls] = map_rows[colInd[j]-startCol];
			target_nzls++;
		}
		tmpRow[target_row+1] = target_nzls;
		target_row++;
			
	}
	int getBelowSize(int *rowInd, int* colInd,int endCol,int nextChareFirstRow, int nextChareLastRow, int* lastRowInd) 
	{
		int size=0;
		for (int i=nextChareFirstRow; i<nextChareLastRow; i++) {
			int j=rowInd[i]+lastRowInd[i];
			while (colInd[j]<endCol) {
				j++;
				size++;
			}
		}
		return size;
	}
	void addBelowRows(int &row_no, int &max_col, int& entries, int* rowInd, int* colInd, double* data, int* tmpRow,
					  int* tmpCol, double* tmpData, bool* tmpDep, int startCol,int endCol, int nextChareFirstRow, int nextChareLastRow, 
					  int* lastRowInd, int* map_rows, vector<chare_deps_str> &chare_deps, row_attr* prev_in_row, int thisChare)
	{
		max_col = 0;
		tmpRow[endCol-startCol+row_no] = entries;
		for (int i=nextChareFirstRow; i<nextChareLastRow; i++) {
			int j=rowInd[i]+lastRowInd[i];
			while (colInd[j]<endCol) {
				tmpData[entries] = data[j];
				tmpCol[entries] = map_rows[colInd[j]-startCol];
				if (tmpCol[entries]>max_col) {
					max_col = tmpCol[entries];
				}
				entries++;
				j++;
			}
			// if row was empty
			if (j==rowInd[i]+lastRowInd[i]) {
				continue;
			}
			if (lastRowInd[i]!=0) {
				tmpDep[endCol-startCol+row_no] = true;
				setDependentChare(chare_deps, prev_in_row[i], thisChare, endCol-startCol+row_no);
			}
			else {
				tmpDep[endCol-startCol+row_no] = false;
			}
			lastRowInd[i] = j-rowInd[i];
			row_attr rtr; rtr.chare= thisChare; rtr.row=row_no;
			prev_in_row[i] = rtr;
			row_no++;
			tmpRow[endCol-startCol+row_no] = entries;
		}
	}
	void setDependentChare(vector<chare_deps_str> &chare_deps, row_attr prev_row, int thisChare, int row_no)
	{
		for (int k=0; k<chare_deps.size(); k++) {
			if (chare_deps[k].chare_no==prev_row.chare) {
				row_attr tmprattr; tmprattr.chare=thisChare; tmprattr.row=row_no;
				chare_deps[k].nextRow[prev_row.row] = tmprattr;
			}
		}
	}
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
};

#include "sparse_solve.def.h"