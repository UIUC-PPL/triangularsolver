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
	
	bool finished;
	bool empty_nondiags;
	
	double *arrived_data;
	bool *arrived_is; // for diag, is item arrived
	
	double* xVal; // value of x
	double* rhs; // right-hand side for each column if diagonal chare
	ArrayMeshStreamer<RowSum, int> * streamer; 
public:
	ColumnsSolve()
	{
		allDone=0;
		belowNumDone=0;
		finished = false;
		empty_nondiags = false;
		first_below_done = false;
		rest_below_done = false;
	}
	ColumnsSolve(CkMigrateMessage *m) {}
	
	void set_input(int num_entries, int num_rows, int num_cols, double m_data[], int m_colInd[], 
				   int m_rowInd[],	bool dep[],	bool m_diag, int m_indep_row_no, int m_first_below_rows,
				   int m_first_below_max_col, int m_rest_below_rows, int m_rest_below_max_col)
	{
		streamer = ((ArrayMeshStreamer<RowSum, int> *)CkLocalBranch(aggregator));
		data = new double[num_entries];
		colInd = new int[num_entries];
		rowInd = new int[num_rows+1];
		row_dep = new bool[num_rows];
		memcpy(data, m_data, num_entries*sizeof(double));
		memcpy(colInd, m_colInd, num_entries*sizeof(int));
		memcpy(rowInd, m_rowInd, (num_rows+1)*sizeof(int));
		memcpy(row_dep, dep, (num_rows)*sizeof(bool));
		nMyCols = num_cols;
		nMyRows = num_rows;
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
		} else {
			nextRow = new row_attr[nMyRows];
		}
		
		arrived_data = new double[nMyRows];
	}
	void set_deps(int dep_size, row_attr deps[])
	{
		memcpy(nextRow, deps, (dep_size)*sizeof(row_attr));
	}

	void set_section(CProxySection_ColumnsSolve nondiags, bool empty_sec)
	{
		lower_section = nondiags;
		empty_nondiags = empty_sec;
	}
	
	void my_indep_compute()
	{
		if (diag && indep_row_no) {
			int i=0;
			// get next chares data faster
			if (first_below_max_col < indep_row_no) {
				for (; i<=first_below_max_col; i++) {
					double val = 0;
					for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
						val += data[j]*xVal[colInd[j]];
					xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
				}
				for (int i=nMyCols; i<nMyCols+first_below_rows; i++) {
					double val=0;
					if (row_dep[i]) {
						if (arrived_is[i]) {
							val = arrived_data[i];
							arrived_is[i] = false;
						}else
							continue;
					}
					
					for (int j=rowInd[i]; j<rowInd[i+1]; j++)
						val += data[j]*xVal[colInd[j]];
					//msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
					streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
					belowNumDone++;
				}
				first_below_done = true;
				//msgPool.flushMsgPool();
			}
			// get other chares data
			if (rest_below_max_col < indep_row_no) {
				for (; i<=rest_below_max_col; i++) {
					double val = 0;
					for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
						val += data[j]*xVal[colInd[j]];
					xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
				}
				for (int i=nMyCols+first_below_rows; i<nMyCols+first_below_rows+rest_below_rows; i++) {
					double val=0;
					if (row_dep[i]) {
						if (arrived_is[i]) {
							val = arrived_data[i];
							arrived_is[i] = false;
						}else 
							continue;
					}
					
					for (int j=rowInd[i]; j<rowInd[i+1]; j++)
						val += data[j]*xVal[colInd[j]];
					//msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
					streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
					belowNumDone++;
				}
				rest_below_done = true;
				//msgPool.flushMsgPool();
			}
			for (; i<indep_row_no; i++) {
				double val = 0;
				for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
					val += data[j]*xVal[colInd[j]];
				xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
			}
			allDone = indep_row_no;
			diag_compute(allDone);
			// send dummy message for sdag while to complete
			//if (finished)
			thisProxy(thisIndex.data[0]).process(RowSum(-1,0));
		}
	}
	void diag_compute(int start)
	{
		// if hadn't started yet
		if (allDone<indep_row_no)
			return;
		int i=start;
		for (; i<=first_below_max_col; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i])
					val = arrived_data[i];
				else
					break;
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
				val += data[j]*xVal[colInd[j]];
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
		}
		if (i> first_below_max_col && !first_below_done) {
			for (int i=nMyCols; i<nMyCols+first_below_rows; i++) {
				double val=0;
				if (row_dep[i]) {
					if (arrived_is[i]) {
						val = arrived_data[i];
						arrived_is[i] = false;
					}else
						continue;
				}
				
				for (int j=rowInd[i]; j<rowInd[i+1]; j++)
					val += data[j]*xVal[colInd[j]];
				// msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
				streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
				belowNumDone++;
			}
			first_below_done = true;
			// msgPool.flushMsgPool();
		}
		for (; i<=rest_below_max_col; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i])
					val = arrived_data[i];
				else
					break;
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
				val += data[j]*xVal[colInd[j]];
			xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
		}
		if (i> rest_below_max_col && !rest_below_done) {
			for (int i=nMyCols+first_below_rows; i<nMyCols+first_below_rows+rest_below_rows; i++) {
				double val=0;
				if (row_dep[i]) {
					if (arrived_is[i]) {
						val = arrived_data[i];
						arrived_is[i] = false;
					}else
						continue;
				}
				
				for (int j=rowInd[i]; j<rowInd[i+1]; j++)
					val += data[j]*xVal[colInd[j]];
				// msgPool.add(nextRow[i-nMyCols].chare, nextRow[i-nMyCols].row, val);
				streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
				belowNumDone++;
			}
			rest_below_done = true;
			// msgPool.flushMsgPool();
		}
		
		for (; i<nMyCols; i++) {
			double val = 0;
			if (row_dep[i-indep_row_no]) {
				if (arrived_is[i])
					val = arrived_data[i];
				else
					break;
			}
			for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
				val += data[j]*xVal[colInd[j]];
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
			if (belowNumDone==first_below_rows+rest_below_rows)
				finished = true;
		}
	}
	void nondiag_compute()
	{
		for (int i=0; i<nMyRows; i++) {
			if (row_dep[i]) {
				continue;
			}
			double val=0;
			for (int j=rowInd[i]; j<rowInd[i+1]; j++)
				val += data[j]*xVal[colInd[j]];
			//msgPool.add(nextRow[i].chare, nextRow[i].row, val);
			streamer->insertData(RowSum(nextRow[i].row, val), nextRow[i].chare);
			allDone++;			
		}
		// msgPool.flushMsgPool();
		if (allDone==nMyRows)
			finished = true;
	}
	void diagReceiveData()
	{
		double val = 0;
		int row = 0;
				if (row==allDone && row<nMyCols) {
					for (int j=rowInd[row]; j<rowInd[row+1]-1; j++)
						val += data[j]*xVal[colInd[j]];
					xVal[row] = (rhs[row]-val)*(data[rowInd[row+1]-1]);
					// first depending row is done
					allDone++;
					diag_compute(allDone);
				}
				// if after diagonal 
				else if (row>=nMyCols && ( (row<nMyCols+first_below_rows && allDone>first_below_max_col)
										 || (row>=nMyCols+first_below_rows && allDone>rest_below_max_col))) {
					for (int j=rowInd[row]; j<rowInd[row+1]; j++)
						val += data[j]*xVal[colInd[j]];
					//msgPool.add(nextRow[row-nMyCols].chare, nextRow[row-nMyCols].row, val);
					streamer->insertData(RowSum(nextRow[row-nMyCols].row, val), nextRow[row-nMyCols].chare);
					belowNumDone++;
				}
				else {
					arrived_data[row] = val;
					arrived_is[row] = true;
				}
//			diag_compute(allDone);
	}
	void nondiagReceiveData(){
		double val = 0;
		int row = 0;
		for (int j=rowInd[row]; j<rowInd[row+1]; j++)
			val += data[j]*xVal[colInd[j]];
		streamer->insertData(RowSum(nextRow[row].row, val), nextRow[row].chare);
		allDone++;
		if (allDone==nMyRows)
			finished = true;
	}
};
