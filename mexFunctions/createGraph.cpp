#include <math.h> /* Needed for the ceil() prototype */
#include "mex.h"
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>

using namespace std;


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// inputs: 
	// index
	// candidateHTrans (nonzeor*(psize*psize)*labels)
	// patchSize
	// m, n
	// structureSzie
	// labels

    if (nrhs != 7) {
		mexErrMsgTxt("7 inputs required: \n index candidateHTrans patchSize structureSize labels \n");
	}

	if (nlhs != 3) {
		mexErrMsgTxt("3 outputs required: \n structure edgePots edgeEnds \n");
	}

    // matlab index
    double *index = (double *) mxGetData(prhs[0]);
    double *candidateHTrans = (double *) mxGetData(prhs[1]);

    int patchSize = mxGetScalar(prhs[2]);
	int m = mxGetScalar(prhs[3]);
	int n = mxGetScalar(prhs[4]);
	int labels = mxGetScalar(prhs[5]);
	int structureSize = mxGetScalar(prhs[6]);

    int half = (patchSize-1)/2;

    cout << "============================================"<<endl;
    cout << "Parameters: "<<endl;
    cout << "patchSize = " << patchSize << endl;
    cout << "m = " << m << " n = " << n <<endl;
    cout << "labels = " << labels <<endl;
    cout << "sturcutreSize = " << structureSize << endl;


    mwIndex *ir,*jc;
    double *sr;
    mwSize nzmax = (mwSize)(structureSize*50);
    plhs[0] = mxCreateSparse(mwSize(structureSize), mwSize(structureSize), nzmax, mxREAL);
    sr  = mxGetPr(plhs[0]);    /* Row indexing      */
    ir = mxGetIr(plhs[0]);	   /* Column count      */
    jc = mxGetJc(plhs[0]);     /* Non-zero elements */

    // first construct stucture
    int nonzero = 0;
    for (int i = 0; i < structureSize; i++) {
    	jc[i] = nonzero;
    	int x0 = index[i]-1;
		int y0 = index[i+structureSize]-1;

		for (int j = i+1; j < structureSize; j++) {
			int xx = index[j]-1;
			int yy = index[j+structureSize]-1;

			if (abs(x0-xx) < patchSize && abs(y0-yy) < patchSize) {
				sr[nonzero] = 1;
				ir[nonzero] = j; //j
				nonzero++;
			}
		}
    }
    jc[structureSize] = nonzero;

    // extract the sparse data
    vector<int> a, b;
    for (int i = 0; i < structureSize; i++) {
    	int stop = jc[i+1];
    	for (int j = jc[i]; j < stop; j++) {
    		a.push_back(ir[j]);     // row
    		b.push_back(i);         // col
    	}
    }

    // c++ index
    // output
	plhs[1] = mxCreateNumericMatrix(nonzero, labels*labels, mxDOUBLE_CLASS, mxREAL);
	double *edgePots = (double *)mxGetData(plhs[1]);

	plhs[2] = mxCreateNumericMatrix(nonzero, 2, mxDOUBLE_CLASS, mxREAL);
	double *edgeEnds = (double *)mxGetData(plhs[2]);


    for (int i = 0 ; i < nonzero; i++) {
    	int tmpMap[m*n];
    	// initialization
    	for (int ii = 0; ii < m*n; ii++) {
    		tmpMap[ii] = 0;
    	}

		int y1 = index[a[i]]-1;
		int x1 = index[a[i]+structureSize]-1;
		int y2 = index[b[i]]-1;
		int x2 = index[b[i]+structureSize]-1;

    	for (int j = -half; j <= half; j++) {
    		for (int k = -half; k <= half; k++) {
    			tmpMap[(x1+k)+(y1+j)*m] ++;
    			tmpMap[(x2+k)+(y2+j)*m] ++;
    		}
    	}

    	// compute the bounding box
    	int min_x = 10000, min_y = 10000, max_x = -10000, max_y = -10000;
    	for (int j = min(y1, y2)-half; j <= max(y1, y2)+half; j++) {
    		for (int k = min(x1, x2)-half; k <= max(x1, x2)+half; k++) {

    			if (tmpMap[k+j*m] == 2) {
    				if (min_x > k) min_x = k;
    				if (min_y > j) min_y = j;
    				if (max_x < k) max_x = k;
    				if (max_y < j) max_y = j;
    			}
    		}
    	}

		
        for (int k1 = 0; k1 < labels; k1++) {
			for (int k2 = 0; k2 < labels; k2++) {

				double sum = 0;
				for (int jy = min_y; jy <= max_y; jy++) {
					for (int jx = min_x; jx <= max_x; jx++) {

						// be careful about matlab index
						double current = candidateHTrans[a[i] + ((jx-x1+half)*patchSize + jy-y1+half)*structureSize + k1*patchSize*patchSize*structureSize];
						double next = candidateHTrans[b[i] + ((jx-x2+half)*patchSize + jy-y2+half)*structureSize + k2*patchSize*patchSize*structureSize];
						double diff = abs(current - next);

						diff *= diff;
						sum += diff;
					}
				}

				// matlab index
				edgePots[i + (k2 + k1*labels)*nonzero] = sum; 
				edgeEnds[i] = a[i]+1; 
				edgeEnds[i+nonzero] = b[i]+1;
    		}
    	}
    }
}