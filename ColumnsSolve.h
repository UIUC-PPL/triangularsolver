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
	ColumnsSolve_SDAG_CODE;
	
	MessagePool msgPool;
	
	int* rowInd;  // index of each row, one more for convenience
	double* data; // data in rows sparse compressed
	int *colInd; // column index of each element
	
	int nMyCols;  // number of columns of this chare
	int nMyRows; // number of rows
	bool diag; // is diagonal Chare
	int indep_row_no; // for diag
	int first_below_rows;
	int first_below_max_col;
	int rest_below_rows;
	int rest_below_max_col;
	bool first_below_done;
	bool rest_below_done;
	
	int allDone; // done rows
	bool* row_dep; // is each row dependent
	CProxySection_ColumnsSolve lower_section; // lower diag chares, for diags
	row_attr* nextRow; // for nondiags
	int belowNumDone;
	
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
public:
	ColumnsSolve():msgPool(thisProxy)
	{
		__sdag_init();
		allDone=0;
		belowNumDone=0;
		xval_got = false;
		finished = false;
		empty_nondiags = false;
		first_below_done = false;
		rest_below_done = false;
	}
	ColumnsSolve(CkMigrateMessage *m):msgPool(thisProxy) {}
	
	void set_input(int num_entries, int num_rows, int num_cols, double m_data[], int m_colInd[], 
				   int m_rowInd[],	bool dep[],	bool m_diag, int m_indep_row_no, int m_first_below_rows,
				   int m_first_below_max_col, int m_rest_below_rows, int m_rest_below_max_col)
	{
		data = m_data;
		nMyCols = num_cols;
		nMyRows = num_rows;
		
		colInd = m_colInd;
		rowInd = m_rowInd;
		row_dep = dep;
		diag = m_diag;
		indep_row_no = m_indep_row_no;
		
		xVal = new double[nMyCols];
		if (diag) {
			rhs = new double[nMyCols];
			for (int i=0; i<nMyCols; i++) {
				rhs[i] = 1.0;
			}
			arrived_is = new bool[nMyRows];
			memset(arrived_is, 0, nMyRows*sizeof(bool));
			first_below_rows = m_first_below_rows;
			first_below_max_col = m_first_below_max_col;
			rest_below_rows = m_rest_below_rows;
			rest_below_max_col = m_rest_below_max_col;
			nextRow = new row_attr[first_below_rows+rest_below_rows];
			arrived_rows = new int[first_below_rows+rest_below_rows];
		} else {
			nextRow = new row_attr[nMyRows];
			arrived_rows = new int[nMyRows];
		}
		
		arrived_data = new double[nMyRows];
		arrived_size = 0;
	}
	// for nondiagonal chares
	void set_deps(DepsMsg* msg)
	{
		int dep_size = diag? first_below_rows+rest_below_rows: nMyRows;
		memcpy(nextRow, msg->deps, (dep_size)*sizeof(row_attr));
		delete msg;
	}
	// for diagonal Chares
	void set_section(CProxySection_ColumnsSolve nondiags, bool empty_sec)
	{
		lower_section = nondiags;
		empty_nondiags = empty_sec;
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
			contribute();
		}
	}
	void start() {
		CkEntryOptions opts;
		opts.setQueueing(CK_QUEUEING_FIFO);
		opts.setPriority(10);
		thisProxy[thisIndex].indep_compute();
	}
	void indep_compute()
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
			contribute();
		}
	}
	void receiveData(DataMsg* msg)
	{
		if (diag) {
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
};
