/*
 author: Ehsan Totoni, totoni2@illinois.edu
 December, 2011
 */


#define MAX(a, b) (a>=b?a:b)
#include <assert.h>
#include <vector>
using namespace std;

/*array [1D]*/
class ColumnsSolve : public CBase_ColumnsSolve
{
	MessagePool msgPool;
	
	int* rowInd;  // index of each row, one more for convenience
	double* data; // data in rows sparse compressed
	int *colInd; // column index of each element
	
	
	int nMyCols;  // number of columns of this chare
	int nMyRows; // number of rows
	// int myColStart; // my first column number global
	// int myRowStart; // my first row number
	bool diag; // is diagonal Chare
	int indep_row_no; // for diag
	int first_below_rows;
	int first_below_max_col;
	int rest_below_rows;
	int rest_below_max_col;
	bool first_below_done;
	bool rest_below_done;
	// bool* below_done;
	
	int allDone; // done rows
	bool* row_dep; // is each row dependent
	CProxySection_ColumnsSolve lower_section; // lower diag chares, for diags
	row_attr* nextRow; // for nondiags
	int belowNumDone;
	
	bool input_got; // input arrived
	bool deps_got; // next row dependencies arrived, for nondiags
	bool section_got; // section proxy arrived, for diags
	bool xval_got; // x values arrived, for nondiagonal chares
	bool finished;
	bool empty_nondiags;
	
	// for nondiags, arrived messages before xVals
	int *arrived_rows;
	double *arrived_data;
	int arrived_size;
	bool *arrived_is; // for diag, is item arrived
	
	double* xVal; // value of x
	double* rhs; // right-hand side for each column if diagonal chare
	int* map; // for checking
public:
	ColumnsSolve():msgPool(thisProxy)
	{
		allDone=0;
		belowNumDone=0;
		input_got = false;
		deps_got = false;
		section_got = false;
		xval_got = false;
		finished = false;
		empty_nondiags = false;
		first_below_done = false;
		rest_below_done = false;
	}
	ColumnsSolve(CkMigrateMessage *m):msgPool(thisProxy) {}
	
	void get_input(InputMsg* msg)
	{
		int data_size = msg->rowInd[msg->num_rows];
		// printf("chare: %d data size: %d\n",thisIndex, data_size);
		data = new double[data_size];
		for (int i=0; i<data_size; i++) {
			data[i] = 1.;
		}
		colInd = new int[data_size];
		rowInd = new int[msg->num_rows+1];
		row_dep = new bool[msg->num_rows];
		// memcpy(data, msg->data, data_size*sizeof(double));
		memcpy(colInd, msg->colInd, data_size*sizeof(int));
		memcpy(rowInd, msg->rowInd, (msg->num_rows+1)*sizeof(int));
		memcpy(row_dep, msg->dep, (msg->num_rows)*sizeof(bool));
		diag = msg->diag;
		nMyCols = msg->num_cols;
		nMyRows = msg->num_rows;
		indep_row_no = msg->indep_row_no;
		xVal = new double[nMyCols];
		if (diag) {
			rhs = new double[nMyCols];
			for (int i=0; i<nMyCols; i++) {
				rhs[i] = i+msg->start_col+1.0;
			}
			arrived_is = new bool[nMyRows];
			memset(arrived_is, 0, nMyRows*sizeof(bool));
			first_below_rows = msg->first_below_rows;
			first_below_max_col = msg->first_below_max_col;
			rest_below_rows = msg->rest_below_rows;
			rest_below_max_col = msg->rest_below_max_col;
			//below_done = new bool[first_below_rows+rest_below_rows];
			//memset(below_done, 0, (first_below_rows+rest_below_rows)*sizeof(bool));
			// analyse_indeps();
			nextRow = new row_attr[first_below_rows+rest_below_rows];
			arrived_rows = new int[first_below_rows+rest_below_rows];
		} else {
			nextRow = new row_attr[nMyRows];
			arrived_rows = new int[nMyRows];
		}
		
		arrived_data = new double[nMyRows];
		arrived_size = 0;
		input_got = true;
		if ((!diag && deps_got) || (diag && section_got && deps_got)) {
			mainProxy.got_input();
		}
		delete msg;
		
	}
	// for nondiagonal chares
	void get_deps(DepsMsg* msg)
	{
		int dep_size = diag? first_below_rows+rest_below_rows: nMyRows;
		memcpy(nextRow, msg->deps, (dep_size)*sizeof(row_attr));
		delete msg;
		deps_got = true;
		if ((!diag &&input_got) || (diag && section_got && input_got)) {
			mainProxy.got_input();
		}
	}
	// for diagonal Chares
	void get_section(CProxySection_ColumnsSolve nondiags, bool empty_sec)
	{
		lower_section = nondiags;
		section_got = true;
		empty_nondiags = empty_sec;
		if (input_got && deps_got) {
			mainProxy.got_input();
		}
	}
	void get_xval(xValMsg* msg)
	{
		memcpy(xVal, msg->xVal, (nMyCols)*sizeof(double));
		delete msg;
		xval_got = true;
		// start computation
		nondiag_compute(0);
		for (int i=0; i<arrived_size; i++) {
			double val = arrived_data[i];
			int row = arrived_rows[i];
			for (int j=rowInd[row]; j<rowInd[row+1]; j++) {
				val += data[j]*xVal[colInd[j]];
			}
			msgPool.add(nextRow[row].chare, nextRow[row].row, val);
			
			allDone++;
		}
		msgPool.flushMsgPool();
		if (allDone==nMyRows && !finished) {
			finished = true;
			// print_check();
			// mainProxy.done();
			contribute();
		}
	}
	void start(DummyMsg* msg)
	{
		if (diag && indep_row_no) {
			int i=0;
			// get next chares data faster
			if (first_below_max_col < indep_row_no) {
				for (; i<=first_below_max_col; i++) {
					double val = 0;
					for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
				}
				for (int i=nMyCols; i<nMyCols+first_below_rows; i++) {
					double val=0;
					if (row_dep[i]) {
						if (arrived_is[i]) {
							val = arrived_data[i];
							arrived_is[i] = false;
						}else {
							continue;
						}
					}
					
					for (int j=rowInd[i]; j<rowInd[i+1]; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
					belowNumDone++;
					// below_done[i-nMyCols] = true;
				}
				first_below_done = true;
				msgPool.flushMsgPool();
			}
			// get other chares data
			if (rest_below_max_col < indep_row_no) {
				for (; i<=rest_below_max_col; i++) {
					double val = 0;
					for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
				}
				for (int i=nMyCols+first_below_rows; i<nMyCols+first_below_rows+rest_below_rows; i++) {
					// if (below_done[i-nMyCols]) {
					//	printf("----??\n");
					//	continue;
					//}
					double val=0;
					if (row_dep[i]) {
						if (arrived_is[i]) {
							val = arrived_data[i];
							arrived_is[i] = false;
						}else {
							continue;
						}
					}
					
					for (int j=rowInd[i]; j<rowInd[i+1]; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
					belowNumDone++;
					// below_done[i-nMyCols] = true;
				}
				rest_below_done = true;
				msgPool.flushMsgPool();
			}
			for (; i<indep_row_no; i++) {
				double val = 0;
				for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
					val += data[j]*xVal[colInd[j]];
				}
				xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
			}
			allDone = indep_row_no;
			diag_compute(allDone);
		}
	}
	void diag_compute(int start)
	{
		// if hadn't started yet
		if (allDone<indep_row_no) {
			return;
		}
		int i=start;
		for (; i<=first_below_max_col; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i]) {
					val = arrived_data[i];
				} else {
					break;
				}
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
				val += data[j]*xVal[colInd[j]];
			}
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
		}
		if (i> first_below_max_col && !first_below_done) {
			for (int i=nMyCols; i<nMyCols+first_below_rows; i++) {
				double val=0;
				if (row_dep[i]) {
					if (arrived_is[i]) {
						val = arrived_data[i];
						arrived_is[i] = false;
					}else {
						continue;
					}
				}
				
				for (int j=rowInd[i]; j<rowInd[i+1]; j++) {
					val += data[j]*xVal[colInd[j]];
				}
				msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
				belowNumDone++;
				// below_done[i-nMyCols] = true;
			}
			first_below_done = true;
			msgPool.flushMsgPool();
		}
		for (; i<=rest_below_max_col; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i]) {
					val = arrived_data[i];
				} else {
					break;
				}
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
				val += data[j]*xVal[colInd[j]];
			}
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
		}
		if (i> rest_below_max_col && !rest_below_done) {
			for (int i=nMyCols+first_below_rows; i<nMyCols+first_below_rows+rest_below_rows; i++) {
				double val=0;
				if (row_dep[i]) {
					if (arrived_is[i]) {
						val = arrived_data[i];
						arrived_is[i] = false;
					}else {
						continue;
					}
				}
				
				for (int j=rowInd[i]; j<rowInd[i+1]; j++) {
					val += data[j]*xVal[colInd[j]];
				}
				msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
				belowNumDone++;
				// below_done[i-nMyCols] = true;
			}
			rest_below_done = true;
			msgPool.flushMsgPool();
		}
		
		for (; i<nMyCols; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i]) {
					val = arrived_data[i];
				} else {
					break;
				}
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
				val += data[j]*xVal[colInd[j]];
			}
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
		}
		allDone = i;
		
		if (allDone==nMyCols && !finished) {
			
			if (!empty_nondiags) {
				empty_nondiags = true; // don't send again!
				// broadcast to other chares of column
				xValMsg* msg = new (nMyCols, 8*sizeof(int)) xValMsg;
				memcpy(msg->xVal, xVal, (nMyCols)*sizeof(double));
				*(int*)CkPriorityPtr(msg) = 0;
				CkSetQueueing(msg, CK_QUEUEING_IFIFO);
				lower_section.get_xval(msg);
			}
			if (belowNumDone==first_below_rows+rest_below_rows) {
				finished = true;	
				// print_check();
				// mainProxy.done();
				contribute();
			}
		}
	}
	void nondiag_compute(int start)
	{
		int i;
		for (int i=start; i<nMyRows; i++) {
			if (row_dep[i]) {
				continue;
			}
			double val=0;
			for (int j=rowInd[i]; j<rowInd[i+1]; j++) {
				val += data[j]*xVal[colInd[j]];
			}
			msgPool.add(nextRow[i].chare, nextRow[i].row, val);
			allDone++;
			
		}
		msgPool.flushMsgPool();
		if (allDone==nMyRows  && !finished) {
			finished = true;
			// print_check();
			// mainProxy.done();
			contribute();
		}
	}
	void receiveData(DataMsg* msg)
	{
		if (diag) {
			// 
			for (int i=0; i<msg->size; i++) {
				double val = msg->data[i];
				int row = msg->rowInd[i];
				if (row==allDone && row<nMyCols) {
					for (int j=rowInd[row]; j<rowInd[row+1]-1; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					xVal[row] = (rhs[row]-val)*(data[rowInd[row+1]-1]);
					// first depending row is done
					allDone++;
					diag_compute(allDone);
				}
				// if after diagonal 
				else if (row>=nMyCols && ( (row<nMyCols+first_below_rows && allDone>first_below_max_col)
										 || (row>=nMyCols+first_below_rows && allDone>rest_below_max_col))) {
					for (int j=rowInd[row]; j<rowInd[row+1]; j++) {
						val += data[j]*xVal[colInd[j]];
					}
					msgPool.add(nextRow[row-nMyCols].chare, nextRow[row-nMyCols].row, val);
					belowNumDone++;
				} 
				else {
					arrived_data[row] = val;
					// if (arrived_is[row]==true) {
					//	printf("index:%d row:%d\n",thisIndex,row);
					//}
					arrived_is[row] = true;
				}	
			}
			msgPool.flushMsgPool();
			diag_compute(allDone);
			
		}
		// nondiag
		else if (xval_got) {
			for (int i=0; i<msg->size; i++) {
				double val = msg->data[i];
				int row = msg->rowInd[i];
				for (int j=rowInd[row]; j<rowInd[row+1]; j++) {
					val += data[j]*xVal[colInd[j]];
				}
				msgPool.add(nextRow[row].chare, nextRow[row].row, val);
				allDone++;
			}
			msgPool.flushMsgPool();
			if (allDone==nMyRows  && !finished) {
				finished = true;
				// print_check();
				// mainProxy.done();
				contribute();
			}
		} // xvalues not arrived yet
		else {
			memcpy(arrived_data+arrived_size, msg->data, (msg->size)*sizeof(double));
			memcpy(arrived_rows+arrived_size, msg->rowInd, (msg->size)*sizeof(int));
			arrived_size += msg->size;
		}
		
		delete msg;
	}
	
	/*	int print_col_done()
	 {
	 for (int i=0; i<nMyCols; i++) {
	 if (i%5==0) {
	 printf("%d:",i);
	 }
	 printf("%d ",doneColumn(i));
	 if (i%30==0 && i) {
	 printf("\n");
	 }
	 }
	 printf("\n");
	 return -1;
	 }
	 */
	void print_check()
	{
		CkPrintf("done: %d\n",thisIndex);
		//if ( thisIndex==44) {
		//	CkPrintf("Chare:%d cols:%d rows:%d done:%d finished:%b\n",thisIndex, nMyCols, nMyRows, finished);
		//}
	}
	// for diags
	void check(xValMsg* msg)
	{
		bool error=false;
		int fails=0;
		for (int i=0; i<nMyCols; i++) {
			if (xVal[map[i]]!=msg->xVal[i]) {
				fails++;
				if (!error) {
					error=true;
					CkPrintf("Chare:%d first error:%d %f %f map:%d\n",thisIndex,i, xVal[map[i]], msg->xVal[i], map[i]);
				}
			}
		}
		if (error) {
			CkPrintf("Chare:%d num errors:%d\n",thisIndex,fails);
		}
		delete msg;
		mainProxy.doneCheck();
	}
	void get_map(MapMsg* msg)
	{
		map = new int[nMyCols];
		memcpy(map, msg->map_rows, (nMyCols)*sizeof(int));
		double *rhst=new double[nMyCols];
		for (int i=0; i<nMyCols; i++) {
			rhst[map[i]] = rhs[i];
		}
		rhs = rhst;
		delete msg;
	}
	void analyse_indeps()
	{
		vector<int> depvec;
		int indeps=0;
		// ignore first
		int firstdep = -1;
		for (int i=0; i<nMyCols; i++) {
			bool dep=false;
			if (row_dep[i]) {
				depvec.push_back(i);
				if (firstdep==-1) {
					firstdep=i;
				}
				continue;
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++) {
				for (int k=0; k<depvec.size(); k++) {
					if (colInd[j]==depvec[k]) {
						dep=true;
						depvec.push_back(i);
						if (firstdep==-1) {
							firstdep=i;
						}
						break;
					}
				}
				if (dep==true) {
					break;
				}
			}
			if (!dep) {
				indeps++;
			}
		}
		CkPrintf("chare:%d indeps:%d of %d rows, first_dep:%d fraction:%f\n",thisIndex, 
				 indeps, nMyCols, firstdep, indeps/(double)nMyCols);
	}
	void status()
	{
		CkPrintf("chare:%d done:%d\n",thisIndex, allDone);
	}
};
