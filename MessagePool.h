/*
 author: Ehsan Totoni, totoni2@illinois.edu
 December, 2011
 */

#define NO_MSG_SEGS 5
#define MSG_SEG_SIZE 100000

class DataMsg : public CMessage_DataMsg
{
public:
	// int *colInd; // column number for each element, or colptr in case of input
    double* data; // data
	int* rowInd;  // row number for each element
	int size;
};

class InputMsg : public CMessage_InputMsg
{
public:
	int* rowInd;  // row ptr
	int *colInd; // column number for each element
	double* data; // data
	bool *dep; // dependency of each column
	int num_rows;
	int start_col;
	int num_cols;
	bool diag; // is diagnoal Chare?
	int indep_row_no; // number of independent rows
	int first_below_rows;
	int first_below_max_col;
	int rest_below_rows;
	int rest_below_max_col;
};

class DummyMsg : public CMessage_DummyMsg
{
};
class DepsMsg : public CMessage_DepsMsg
{
public:
	row_attr* deps; // solution x
};

class xValMsg : public CMessage_xValMsg
{
public:
	double *xVal; // solution x
};

class MessagePool {
	CProxy_ColumnsSolve thisProxy;
	
	
	// message pool data structures
	double* out_data; // store data to be sent out
	// pre-allocated segments of size MSG_SEG_SIZE for NO_MSG_SEGS PEs
	int* out_rowInd; // row index of each element
	int *out_colsInd; // column index of each element in data array
	int* out_size;    // data size of each processor
	int* index_out; // mapping of processor number to out data structure's index
	int lastChare; // last chare accessed, for faster search
	int lastChare_ind; // mapped index
	int no_pe_segs; // number of segments occupied

public:
	MessagePool(CProxy_ColumnsSolve p)
	{
		thisProxy = p;
		out_data = new double[NO_MSG_SEGS*MSG_SEG_SIZE];
		out_rowInd = new int[NO_MSG_SEGS*MSG_SEG_SIZE];
//		out_colsInd = new int[NO_MSG_SEGS*MSG_SEG_SIZE];
		out_size = new int[NO_MSG_SEGS];
		memset(out_size, 0, NO_MSG_SEGS*sizeof(int));
		index_out = new int[NO_MSG_SEGS];
		for (int i=0; i<NO_MSG_SEGS; i++) {
			index_out[i] = -1;
		}
		lastChare = -1;
		lastChare_ind = -1;
		no_pe_segs = 0;
	}
	
	void add(int chare_no, int row, double val) {
		
		int chare_index;
		int out_ind;
		// find PE
		if (chare_no==lastChare && no_pe_segs) { // last chare and not flushed
			chare_index = lastChare_ind;
		} else {
			int i;
			for (i=0; i<no_pe_segs; i++) {
				if (index_out[i]==chare_no) {
					chare_index = i;
					break;
				}
			}
			// if PE is not allocated
			if (i==no_pe_segs) {
				chare_index = allocMsgSeg();
				index_out[chare_index] = chare_no;
			}
		}
		lastChare = chare_no;
		lastChare_ind = chare_index;
		out_ind = out_size[chare_index]++;
		out_ind += chare_index*MSG_SEG_SIZE;
		out_data[out_ind] = val;
		out_rowInd[out_ind] = row;
		// out_colsInd[out_ind] = col;//convertToGlobalColumn(col);
		// check for full buffer
		if (out_size[lastChare_ind]==MSG_SEG_SIZE) {
			flushMsgPool();
		}
	}
	int allocMsgSeg()
	{
		int PE_index = -1;
		if (no_pe_segs!=NO_MSG_SEGS) {
			/// initialize---------
			PE_index = no_pe_segs++;
		} else {
			flushMsgPool();
			PE_index = 0;
			no_pe_segs=1;
		}
		out_size[PE_index] = 0;
		return PE_index;
	}
	void flushMsgPool() 
	{
		for (int i=0; i<no_pe_segs; i++) {
			// one integer for priority
			DataMsg *msg = new ( out_size[i], out_size[i], 8*sizeof(int)) DataMsg;
			// DataMsg *msg = new ( out_size[i], out_size[i]) DataMsg;
			msg->size = out_size[i];
			memcpy(msg->data, &out_data[i*MSG_SEG_SIZE], out_size[i]*sizeof(double));
			memcpy(msg->rowInd, &out_rowInd[i*MSG_SEG_SIZE],out_size[i]*sizeof(int));
			// assuming first has smallest row index, which is not general!!!
			// *(int*)CkPriorityPtr(msg) = msg->rowInd[0];
			//*(int*)CkPriorityPtr(msg) = ~index_out[i];
//			*(int*)CkPriorityPtr(msg) = ~msg->colInd[msg->size-1];
			*(int*)CkPriorityPtr(msg) = 0;
			CkSetQueueing(msg, CK_QUEUEING_IFIFO);
			thisProxy[index_out[i]].receiveData(msg);
		}
		no_pe_segs=0;
	}

};
