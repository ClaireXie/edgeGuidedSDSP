#include <math.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Variables */
	int n, s, s1, s2, n1, n2, e,y1,y2,
			nNodes, nEdges, maxState,
			*y,
			*edgeEnds;
	
	double *nodePot, *edgePot, *modifiedNP, *modifiedEP, *edgeWeights,tmp;
	
	/* Input */
	
	y = (int*)mxGetPr(prhs[0]);
	s1 = (int)mxGetScalar(prhs[1]);
	s2 = (int)mxGetScalar(prhs[2]);
	nodePot = mxGetPr(prhs[3]);
	edgePot = mxGetPr(prhs[4]);
	edgeEnds = (int*)mxGetPr(prhs[5]);
	modifiedNP = mxGetPr(prhs[6]);
	modifiedEP = mxGetPr(prhs[7]);
	
	if (!mxIsClass(prhs[0],"int32")||!mxIsClass(prhs[1],"int32")||!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[5],"int32"))
		mexErrMsgTxt("y, s1, s2, edgeEnds must be int32");
	
	/* Compute Sizes */
	
	nNodes = mxGetDimensions(prhs[3])[0];
	maxState = mxGetDimensions(prhs[3])[1];
	nEdges = mxGetDimensions(prhs[5])[0];
	
	/* Make conditional node energies */
	for(n=0; n < nNodes; n++) {
		if (y[n]-1==s1) {
			modifiedNP[n] = nodePot[n+nNodes*s1];
			modifiedNP[n+nNodes] = nodePot[n+nNodes*s2];
		}
		else {
			modifiedNP[n] = nodePot[n+nNodes*s1];
			modifiedNP[n+nNodes] = nodePot[n+nNodes*(y[n]-1)];
		}
	}
	for(e=0; e < nEdges; e++) {
		n1 = edgeEnds[e]-1;
		n2 = edgeEnds[e+nEdges]-1;
		y1 = y[n1]-1;
		y2 = y[n2]-1;
		if(y1==s1)
			y1 = s2;
		if(y2==s1)
			y2 = s2;
		
		/* Make conditional edge energies */
		modifiedEP[4*e] = -log(edgePot[s1+maxState*(s1+maxState*e)]); /* (1,1) */
		modifiedEP[1+4*e] = -log(edgePot[y1+maxState*(s1+maxState*e)]); /* (2,1) */
		modifiedEP[2+4*e] = -log(edgePot[s1+maxState*(y2+maxState*e)]); /* (1,2) */
		modifiedEP[3+4*e] = -log(edgePot[y1+maxState*(y2+maxState*e)]); /* (2,2) */
		
		/* Modify to be sub-modular */
		if (y[n1]-1 == s1 && y[n2]-1 == s2) {
			tmp = modifiedEP[2+4*e]+modifiedEP[1+4*e]-modifiedEP[3+4*e];
			if (modifiedEP[4*e] > tmp)
				modifiedEP[4*e] = tmp;
		}
		else if (y[n1]-1 == s1) {
			tmp = modifiedEP[4*e]+modifiedEP[3+4*e]-modifiedEP[2+4*e];
			if (modifiedEP[1+4*e] < tmp)
				modifiedEP[1+4*e] = tmp;
		}
		else if (y[n2]-1 == s1) {
			tmp = modifiedEP[4*e]+modifiedEP[3+4*e]-modifiedEP[1+4*e];
			if (modifiedEP[2+4*e] < tmp)
				modifiedEP[2+4*e] = tmp;
		}
		else {
			tmp = modifiedEP[2+4*e]+modifiedEP[1+4*e]-modifiedEP[4*e];
			if (modifiedEP[3+4*e] > tmp)
				modifiedEP[3+4*e] = tmp;
		}
		modifiedEP[4*e] = exp(-modifiedEP[4*e]);
		modifiedEP[1+4*e] = exp(-modifiedEP[1+4*e]);
		modifiedEP[2+4*e] = exp(-modifiedEP[2+4*e]);
		modifiedEP[3+4*e] = exp(-modifiedEP[3+4*e]);
	}
}
