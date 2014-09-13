#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Variables */
    
    int n, s, f, n1, n2, s1, s2, i,e,
            nInstances, nNodes, nNodeFeatures, nEdges, maxState, nEdgeFeatures,
            *edgeEnds, *nStates, *nodeMap, *edgeMap, *Y;
    double obs, *g, *nodeBel, *edgeBel, *Xnode, *Xedge;
    
    /* Input */
    g = mxGetPr(prhs[0]);
    i = (int)mxGetScalar(prhs[1])-1;
    nodeBel = mxGetPr(prhs[2]);
    edgeBel = mxGetPr(prhs[3]);
    edgeEnds = (int*)mxGetPr(prhs[4]);
    nStates = (int*)mxGetPr(prhs[5]);
    nodeMap = (int*)mxGetPr(prhs[6]);
    edgeMap = (int*)mxGetPr(prhs[7]);
    Xnode = mxGetPr(prhs[8]);
    Xedge = mxGetPr(prhs[9]);
    Y = (int*)mxGetPr(prhs[10]);
	
	if (!mxIsClass(prhs[1],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32")||!mxIsClass(prhs[7],"int32")||!mxIsClass(prhs[10],"int32"))
		mexErrMsgTxt("edgeEnds, nStates, nodeMap, edgeMap, i, Y must be int32");
    
    /* Compute Sizes */
    nInstances = mxGetDimensions(prhs[10])[0];
    nNodeFeatures = mxGetDimensions(prhs[8])[1];
    nNodes = mxGetDimensions(prhs[2])[0];
    nEdgeFeatures = mxGetDimensions(prhs[9])[1];
    nEdges = mxGetDimensions(prhs[4])[0];
    maxState = getMaxState(nStates, nNodes);
    
    for(n = 0; n < nNodes; n++) {
        for(s = 0; s < nStates[n]; s++) {
            for(f = 0; f < nNodeFeatures; f++) {
                if(nodeMap[n + nNodes*(s + maxState*f)] > 0) {
                    if(s == Y[i + nInstances*n]-1) {
                        obs = 1;
                    }
                    else {
                        obs = 0;
                    }
                    g[nodeMap[n + nNodes*(s + maxState*f)]-1] += Xnode[i + nInstances*(f + nNodeFeatures*n)]*(nodeBel[n + nNodes*s] - obs);
                }
            }
        }
    }
     
   
    for(e = 0; e < nEdges; e++) {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        for (s1 = 0; s1 < nStates[n1]; s1++) {
            for (s2 = 0; s2 < nStates[n2]; s2++) {
                for (f = 0; f < nEdgeFeatures; f++) {
                    if (edgeMap[s1 + maxState*(s2 + maxState*(e + nEdges*f))] > 0) {
                        if (s1 == Y[i + nInstances*n1]-1 && s2 == Y[i + nInstances*n2]-1) {
                            obs = 1;
                        }
                        else {
                            obs = 0;
                        }
                        g[edgeMap[s1 + maxState*(s2 + maxState*(e + nEdges*f))]-1] += Xedge[i + nInstances*(f + nEdgeFeatures*e)]*(edgeBel[s1 + maxState*(s2 + maxState*e)] - obs);
                    }
                }
            }
        }
    }
    
}
