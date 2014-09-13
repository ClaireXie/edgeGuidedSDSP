#include <math.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Variables */
	int n,n1,n2,e,nNodes,maxState,nEdges,*edgeEnds,s,*y,*swapPositions,dims[3];
	double *nodePot,*edgePot,*newNodePot,*newEdgePot;
	
	/* Inputs */
	nodePot = mxGetPr(prhs[0]);
	edgePot = mxGetPr(prhs[1]);
	edgeEnds = (int*)mxGetPr(prhs[2]);
	s = (int)mxGetScalar(prhs[3]);
	y = (int*)mxGetPr(prhs[4]);
	swapPositions = (int*)mxGetPr(prhs[5]);

	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32"))
		mexErrMsgTxt("edgeEnds, s, y, swapPositions must be int32");
	
	nNodes = mxGetDimensions(prhs[0])[0];
	maxState = mxGetDimensions(prhs[0])[1];
	nEdges = mxGetDimensions(prhs[2])[0];

	/* Outpus */
	plhs[0] = mxCreateDoubleMatrix(nNodes,2,mxREAL);
	newNodePot = mxGetPr(plhs[0]);
	dims[0] = 2;
	dims[1] = 2;
	dims[2] = nEdges;
	plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	newEdgePot = mxGetPr(plhs[1]);
	
	for(n=0;n<nNodes;n++) {
		newNodePot[n] = nodePot[n+nNodes*s];
		newNodePot[n+nNodes] = nodePot[n+nNodes*(y[swapPositions[n]-1]-1)];
	}
	
	for(e=0;e<nEdges;e++) {
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		newEdgePot[0+2*(0+2*e)] = edgePot[s+maxState*(s+maxState*e)];
		newEdgePot[0+2*(1+2*e)] = edgePot[s+maxState*(y[swapPositions[n2]-1]-1+maxState*e)];
		newEdgePot[1+2*(0+2*e)] = edgePot[y[swapPositions[n1]-1]-1+maxState*(s+maxState*e)];
		newEdgePot[1+2*(1+2*e)] = edgePot[y[swapPositions[n1]-1]-1+maxState*(y[swapPositions[n2]-1]-1+maxState*e)];
	}
}
