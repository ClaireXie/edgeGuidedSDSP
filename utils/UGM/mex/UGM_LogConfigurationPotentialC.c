#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Variables */
    
    int n, n1, n2, e,
            nNodes, nEdges, maxState,
            *edgeEnds, *y;
    double *nodePot, *edgePot, *logPot;
    
    /* Input */
    y = (int*)mxGetPr(prhs[0]);
    nodePot = mxGetPr(prhs[1]);
    edgePot = mxGetPr(prhs[2]);
    edgeEnds = (int*)mxGetPr(prhs[3]);
	
	if (!mxIsClass(prhs[0],"int32")||!mxIsClass(prhs[3],"int32"))
		mexErrMsgTxt("y, edgeEnds must be int32");
    
    /* Compute Sizes */
    nNodes = mxGetDimensions(prhs[1])[0];
    maxState = mxGetDimensions(prhs[1])[1];
    nEdges = mxGetDimensions(prhs[3])[0];
    
    /* Output */
	plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
    logPot = mxGetPr(plhs[0]);
    
    *logPot = 0;
    
    for(n = 0; n < nNodes; n++) {
        *logPot += log(nodePot[n + nNodes*(y[n]-1)]);
    }
    for(e = 0; e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        *logPot += log(edgePot[y[n1]-1 + maxState*(y[n2]-1 + maxState*e)]);
    }
    
    
}
