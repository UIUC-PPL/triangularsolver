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

extern void readInput(char * fileName,double * &val, int * &rowind,int * &colptr,int & m, int&  n, int & nzl);

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
	bool print_seq;
public:
	Main(CkArgMsg* m)
	{
		check=false;
		print_seq = false;
		//Process command-line arguments
		nElements=2;
		if(m->argc >1 ) nElements=atoi(m->argv[1]);
		if(m->argc >2 ) strcpy(fileName, m->argv[2]);
		if (m->argc>3) {
			if(!strcmp(m->argv[3],"-seq")) {
				if (m->argc>4 && !strcmp(m->argv[4],"-p")) {
					print_seq = true;
				}
				seq_solve();
				CkExit();
			}
			if (!strcmp(m->argv[3],"-check")) {
				seq_solve();
				check = true;
			}
		}
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
		CkPrintf("Running ColumnsSolve on %d processors for %d elements\n",
				 CkNumPes(),nElements);
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
		// arr.status();

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

			// if there are nondiagonal allocate dependencies
			// if (firstBelowChunkRows!=0 || restBelowChunkRows!=0) {
				chare_deps_str new_chare; new_chare.size=firstBelowChunkRows+restBelowChunkRows; new_chare.chare_no=i;
				new_chare.nextRow=new row_attr[firstBelowChunkRows+restBelowChunkRows];
				chare_deps.push_back(new_chare);
			// }
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
			if (check) {
				MapMsg* mmsg = new (endCol-startCol) MapMsg;
				memcpy(mmsg->map_rows, map_rows, (endCol-startCol)*sizeof(int));
				arr[i].get_map(mmsg);
			}
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
						// tmpData[entries] = data[j];
						tmpCol[entries] = map_rows[colsInd[j]-startCol];
						// if (tmpCol[entries]>restBelowChunkMaxCol) {
						//	restBelowChunkMaxCol = tmpCol[entries];
						//}
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
	void done(void)
	{
		static int count=0;
		count++;
		if(count==total_chares) {
			CkPrintf("All done in %f\n",CmiWallTimer()-start_time);
			if (check) {
				int bsize=total_columns/nElements;
				for (int i=0; i<nElements-1; i++) {
					xValMsg* msg = new (bsize) xValMsg;
					memcpy(msg->xVal, &xVal[bsize*i], bsize*sizeof(double));
					arr[i].check(msg);
				}
				int before = bsize*(nElements-1);
				bsize = total_columns-before;
				xValMsg* msg = new (bsize) xValMsg;
				memcpy(msg->xVal, &xVal[before], bsize*sizeof(double));
				arr[nElements-1].check(msg);
				CkPrintf("checking solution:\n");
				check=false;
				count=0;
			} else {
				CkExit();
			}
		}
	}
	void reportIn()
	{
		CkPrintf("All done in %f\n",CmiWallTimer()-start_time);
		if (check) {
			int bsize=total_columns/nElements;
			for (int i=0; i<nElements-1; i++) {
				xValMsg* msg = new (bsize) xValMsg;
				memcpy(msg->xVal, &xVal[bsize*i], bsize*sizeof(double));
				arr[i].check(msg);
			}
			int before = bsize*(nElements-1);
			bsize = total_columns-before;
			xValMsg* msg = new (bsize) xValMsg;
			memcpy(msg->xVal, &xVal[before], bsize*sizeof(double));
			arr[nElements-1].check(msg);
			CkPrintf("checking solution:\n");
			check=false;
		} else {
			CkExit();
		}
		
	}
	void doneCheck(void)
	{
		static int count=0;
		count++;
		if(count==nElements) {
			CkExit();
		}
	}
	void got_input(void)
	{
		static int count=0;
		count++;
		if(count==total_chares) {
			//CkEntryOptions opts;
			//opts.setPriority(1);
			start_time = CmiWallTimer();
			for(int i=0; i<nElements; i++)
			{
			DummyMsg *msg = new (8*sizeof(int)) DummyMsg;
			*(int*)CkPriorityPtr(msg) = 10000;
			CkSetQueueing(msg, CK_QUEUEING_IFIFO);
				arr[i].start(msg);
			}
		}
	};
	void seq_solve()
	{
		double* data; // data in row sparse compressed
		int* rowInd;  // row index of each element
		int *colsInd; // index of each column in data array, one more than columns for convenience
		// read input
		int m,n, nzl;
		readInput(fileName, data, rowInd, colsInd ,m, n, nzl);
		double *rhs = new double[m];
		for (int i=0; i<m; i++) {
			rhs[i]=i+1.0;
		}
		xVal = new double[n];
		memset(xVal, 0, n*sizeof(double));
		
		double start=CmiWallTimer();
		
		for (int i=0; i<n; i++) {
			double val=0;
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
				val += xVal[colsInd[j]]*data[j];
			}
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
			if (print_seq)
				printf("x[%d]=%f ",i,xVal[i]);
		}
		if (print_seq)
			printf("\n");
		printf("seq time: %f\n",CmiWallTimer()-start);
		delete[] data;
		delete[] rowInd;
		delete[] colsInd;
		delete[] rhs;
	}
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
};


#include "sparse_solve.def.h"
