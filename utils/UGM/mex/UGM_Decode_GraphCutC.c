#include <math.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Variables */
	int n1,n2,e,nNodes,nEdges,*edgeEnds;
	double *nodeEnergy,*edgeEnergy;
	
	/* Inputs */
	nodeEnergy = mxGetPr(prhs[0]);
	edgeEnergy = mxGetPr(prhs[1]);
	edgeEnds = (int*)mxGetPr(prhs[2]);
	
	if (!mxIsClass(prhs[2],"int32"))
        mexErrMsgTxt("edgeEnds must be int32");
	
	nNodes = mxGetDimensions(prhs[0])[0];
	nEdges = mxGetDimensions(prhs[2])[0];
	
	for(e=0;e<nEdges;e++) {
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		nodeEnergy[n1+nNodes] += edgeEnergy[1 + 2*(0 + 2*e)] - edgeEnergy[0 + 2*(0 + 2*e)];
		nodeEnergy[n2+nNodes] += edgeEnergy[1 + 2*(1 + 2*e)] - edgeEnergy[1 + 2*(0 + 2*e)]; 
	}
}
