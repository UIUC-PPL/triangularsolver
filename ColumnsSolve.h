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
	
	int nMyCols, nMyRows;  // number of columns and rows of this chare
	bool diag; // is diagonal Chare
	int indep_row_no, first_below_rows, first_below_max_col, rest_below_rows, rest_below_max_col;
	bool first_below_done, rest_below_done;
	
	int allDone, belowNumDone; // done rows
	bool* row_dep; // is each row dependent
	CProxySection_ColumnsSolve lower_section; // lower diag chares, for diags
	row_attr* nextRow; // for nondiags
	
	bool finished, empty_nondiags, xval_got;
	xValMsg* xval_msg;
	
	int *arrived_rows;
	double *arrived_data;
	int arrived_size;
	bool *arrived_is; // for diag, is item arrived
	
	double* xVal; // value of x
	double* rhs; // right-hand side for each column if diagonal chare
	ArrayMeshStreamer<RowSum, int> * streamer; 
public:
	ColumnsSolve()
	{
		arrived_size = 0;
		allDone=0;
		belowNumDone=0;
		finished = false;
		empty_nondiags = false;
		first_below_done = false;
		rest_below_done = false;
		xval_got = false;
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
			arrived_rows = new int[first_below_rows+rest_below_rows];
			if (!empty_nondiags) {
				xval_msg = new (nMyCols, 8*sizeof(int)) xValMsg;
				*(int*)CkPriorityPtr(xval_msg) = 0;
				CkSetQueueing(xval_msg, CK_QUEUEING_IFIFO);
				delete[] xVal;
				xVal = xval_msg->xVal;
			}
			
		} else {
			nextRow = new row_attr[nMyRows];
			arrived_rows = new int[nMyRows];
		}
		
		arrived_data = new double[nMyRows];
	}
	void set_deps(int dep_size, row_attr deps[])
	{
		memcpy(nextRow, deps, (dep_size)*sizeof(row_attr));
	}

	void set_section(CProxySection_ColumnsSolve nondiags, bool empty_sec, CkGroupID mCastGrpId)
	{
		lower_section = nondiags;
		empty_nondiags = empty_sec;
		if (!empty_nondiags)
			lower_section.ckSectionDelegate(CProxy_CkMulticastMgr(mCastGrpId).ckLocalBranch());
	}
	
	void indep_compute(int a)
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
					streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
					belowNumDone++;
				}
				first_below_done = true;
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
					streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
					belowNumDone++;
				}
				rest_below_done = true;
			}
			for (; i<indep_row_no; i++) {
				double val = 0;
				for (int j=rowInd[i]; j<rowInd[i+1]-1; j++)
					val += data[j]*xVal[colInd[j]];
				xVal[i] = (rhs[i]-val)*(data[rowInd[i+1]-1]);
			}
			allDone = indep_row_no;
			diag_compute(allDone);
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
					} else
						continue;
				}
				
				for (int j=rowInd[i]; j<rowInd[i+1]; j++)
					val += data[j]*xVal[colInd[j]];
				streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
				belowNumDone++;
			}
			first_below_done = true;
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
				streamer->insertData(RowSum(nextRow[i-nMyCols].row, val), nextRow[i-nMyCols].chare);
				belowNumDone++;
			}
			rest_below_done = true;
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
				lower_section.get_xval(xval_msg);
			}
			if (belowNumDone==first_below_rows+rest_below_rows)
				finish();
		}
	}
	void get_xval(xValMsg* msg)
	{
//		CmiReference(UsrToEnv(msg));
		xval_msg = msg;
		xVal = msg->xVal;
		xval_got = true;
		// start computation
		nondiag_compute();
		for (int i=0; i<arrived_size; i++) {
			double val = arrived_data[i];
			int row = arrived_rows[i];
			for (int j=rowInd[row]; j<rowInd[row+1]; j++)
				val += data[j]*xVal[colInd[j]];
			streamer->insertData(RowSum(nextRow[row].row, val), nextRow[row].chare);
			allDone++;
		}
		if (allDone==nMyRows && !finished)
			finish();
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
			streamer->insertData(RowSum(nextRow[i].row, val), nextRow[i].chare);
			allDone++;			
		}
		if (allDone==nMyRows && !finished)
			finish();
	}
	inline virtual void process(const RowSum &rs){
		double val = rs.val;
		int row = rs.row;
		if (diag) {
				if (row==allDone && row<nMyCols) {
					for (int j=rowInd[row]; j<rowInd[row+1]-1; j++)
						val += data[j]*xVal[colInd[j]];
					xVal[row] = (rhs[row]-val)*(data[rowInd[row+1]-1]);
					// first depending row is done
					allDone++;
				}
				// if after diagonal 
				else if (row>=nMyCols && ( (row<nMyCols+first_below_rows && allDone>first_below_max_col)
										  || (row>=nMyCols+first_below_rows && allDone>rest_below_max_col))) {
					for (int j=rowInd[row]; j<rowInd[row+1]; j++)
						val += data[j]*xVal[colInd[j]];
					streamer->insertData(RowSum(nextRow[row-nMyCols].row, val), nextRow[row-nMyCols].chare);
					belowNumDone++;
				} else {
					arrived_data[row] = val;
					arrived_is[row] = true;
				}	
				diag_compute(allDone);
			}
		// nondiag
		else if (xval_got) {
				for (int j=rowInd[row]; j<rowInd[row+1]; j++)
					val += data[j]*xVal[colInd[j]];
				streamer->insertData(RowSum(nextRow[row].row, val), nextRow[row].chare);
				allDone++;
			if (allDone==nMyRows && !finished)
				finish();
		} else {
			arrived_data[arrived_size]=val;
			arrived_rows[arrived_size++]=row;
		}
	}
	void finish() {
//		if (!diag)
//			CmiFree(UsrToEnv(xval_msg));
		finished = true;
		streamer->done();
	}
};
