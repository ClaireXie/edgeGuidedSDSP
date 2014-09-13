#include <math.h>
#include "mex.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	/* Variables */
    int n,n1,n2,s,s1,s2,nodeNum,e,Vind,edgeNum,
            nNodes,nFree,nEdges,maxState,dims[2],dims3[3],nInducedEdges,
            *nStates,*edgeEnds,*V,*E,*clamped,*newEdgeEnds,*newNstates,*nodeMap,*edgeMap;
    double *nodePot, *edgePot,*newNodePot,*newEdgePot;
	
	/* Inputs */
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    nStates = (int*)mxGetPr(prhs[2]);
    edgeEnds = (int*)mxGetPr(prhs[3]);
    V = (int*)mxGetPr(prhs[4]);
    E = (int*)mxGetPr(prhs[5]);
    clamped = (int*)mxGetPr(prhs[6]);
	
	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32"))
        mexErrMsgTxt("edgeEnds, nStates, V, E, clamped must be int32");
    
    /* Sizes */
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[3])[0];
    
    /* Compute size of induced sub-graph */
    nFree = 0;
    for(n=0;n < nNodes;n++) {
        if(clamped[n]==0)
            nFree++;
    }
    nInducedEdges = 0;
    for(e=0;e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        if(clamped[n1]==0 && clamped[n2]==0)
            nInducedEdges++;
    }
            
    /* Output */
    plhs[0] = mxCreateDoubleMatrix(nFree,maxState,mxREAL);
	dims[0] = nNodes;
	dims[1] = 1;
    plhs[1] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
	dims[0] = nFree;
    plhs[2] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
	dims[0] = maxState;
    dims[1] = maxState;
    dims[2] = nInducedEdges;
    plhs[3] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
	dims[0] = nInducedEdges;
	dims[1] = 2;
	plhs[4] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
	dims[0] = nEdges;
	dims[1] = 1;
	plhs[5] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    newNodePot = mxGetPr(plhs[0]);
    nodeMap = mxGetData(plhs[1]);
    newNstates = mxGetData(plhs[2]);
    newEdgePot = mxGetPr(plhs[3]);
	newEdgeEnds = mxGetData(plhs[4]);
    edgeMap = mxGetData(plhs[5]);
	
    nodeNum = 0;
    for(n=0;n < nNodes;n++) {
        if(clamped[n]==0) {
            
            /* Grab node potential */
            for(s=0;s < nStates[n];s++) {
                newNodePot[nodeNum + nFree*s] = nodePot[n + nNodes*s];
            }
            
            /* Absorb edge potentials from clamped neighbors */
            for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++) {
                e = E[Vind]-1;
                n1 = edgeEnds[e]-1;
                n2 = edgeEnds[e+nEdges]-1;
                
                if(n==n1) {
                    if(clamped[n2]!=0) {
                        for(s=0;s<nStates[n];s++) {
                            newNodePot[nodeNum + nFree*s] *= edgePot[s+maxState*(clamped[n2]-1 + maxState*e)];
                        }
                    }
                }
                else {
                    if(clamped[n1]!=0) {
                        for(s=0;s<nStates[n];s++) {
                            newNodePot[nodeNum + nFree*s] *= edgePot[clamped[n1]-1+maxState*(s + maxState*e)];
                        }
                    }
                }
            }
            nodeMap[n] = nodeNum+1;
            newNstates[nodeNum] = nStates[n];
            nodeNum++;
        }
    }
    
    /* Grab edges in subgraph */
    edgeNum = 0;
    for(e=0;e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        if(clamped[n1]==0 && clamped[n2]==0) {
			for(s1=0;s1<nStates[n1];s1++) {
                for(s2=0;s2<nStates[n2];s2++) {
                    newEdgePot[s1+maxState*(s2+maxState*edgeNum)] = edgePot[s1+maxState*(s2+maxState*e)];
                }
            }
			newEdgeEnds[edgeNum] = nodeMap[n1];
			newEdgeEnds[edgeNum+nInducedEdges] = nodeMap[n2];
			edgeMap[e] = edgeNum+1;
            edgeNum++;
        }
    }
}
