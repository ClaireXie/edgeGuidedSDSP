#include "ANN.h"
#include <vector>
#include <fstream>
#include <iostream>
#include "mex.h"

using namespace std;

/* Retrieve data from mexArray
 * @params ANNpointArray data structure (num*dim)
 * @params pointer to the mexArray data
 * @params number of points
 * @params dimension of the data
 */ 
void retrieve_data(ANNpointArray &dataPts, double *data, int num, int dim) {
	for (int i = 0; i<num; i++) {
		for (int j = 0; j<dim; j++) {
			dataPts[i][j] = data[i + j*num];
		}
	}
}


void mexFunction(int nlhs, mxArray * plhs[], int nrhs, const mxArray * prhs[]){
	
	/* data interface
	* @params training data (N1*DIM)
	* @params query data (not limited to a single point) (N2*DIM)
	* @params number of candidates
	* @params verbose (optional)
	*/
	if (nrhs < 3) {
		mexErrMsgTxt("3 inputs required: data_train (N1*dim), data_query(N2*dim), k, (verbose)");
	}

	bool verbose = false;
	if (nrhs == 4) {
		verbose = mxGetScalar(prhs[3]);
	}

	if (!(mxIsDouble(prhs[0])) || !(mxIsDouble(prhs[1]))) {
		mexErrMsgTxt("Input data should be double.");
	}

	if (mxGetN(prhs[0]) != mxGetN(prhs[1])) {
		mexErrMsgTxt("Dimension of training data and query data should be equal.");
	}

	double* train_data = mxGetPr(prhs[0]);
	double* query_data = mxGetPr(prhs[1]);

	int dim = mxGetN(prhs[0]);
	int num = mxGetM(prhs[0]);
	int num_query = mxGetM(prhs[1]);

	if (verbose) {
		cout << "Load " << num << " Points. " << endl;
		cout << "Load " << num_query << " Query points. " << endl;
	}

	int k = mxGetScalar(prhs[2]);

	// retrieve the data to ann data type
	ANNpointArray dataPts;
	ANNpointArray queryPt;
	vector<ANNidxArray>	nnIdx;	
	vector<ANNdistArray> dists;	
	ANNkd_tree*	kdTree;

	dataPts = annAllocPts(num, dim);	
	queryPt = annAllocPts(num_query, dim);	

	retrieve_data(dataPts, train_data, num, dim);
	retrieve_data(queryPt, query_data, num_query, dim);

	// build the kd tree
	if (verbose) cout << "Building the kd tree..." << endl;
	kdTree = new ANNkd_tree(dataPts, num, dim);	

	// knn query
	if (verbose) cout << "knn searching..." << endl;
	for (int i = 0; i<num_query; i++) {
		ANNidxArray	nnIdx0 = new ANNidx[k];	
		ANNdistArray dists0 = new ANNdist[k];
		kdTree->annkSearch(	queryPt[i], k, nnIdx0, dists0, 0);				

		nnIdx.push_back(nnIdx0);
		dists.push_back(dists0);
	}

	// dump the data to matlab
	plhs[0] = mxCreateDoubleMatrix(nnIdx.size(), k, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(dists.size(), k, mxREAL);
    double* indexes = mxGetPr(plhs[0]);
    double* distances   = mxGetPr(plhs[1]);

    for (int i=0; i < nnIdx.size(); i++){
    	for (int j=0; j < k; j++){
    		indexes[ i+ nnIdx.size()*j] = nnIdx[i][j] + 1;
    		distances[ i+ dists.size()*j] = dists[i][j];
    	}
    }  
}
