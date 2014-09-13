#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* Variables */
    
    int n, s, f, n1, n2, s1, s2, i,e,
            nInstances, nNodes, nNodeFeatures, nEdges, maxState, nEdgeFeatures,
            *edgeEnds, *nStates, *nodeMap, *edgeMap, *Y, *V, *E, Vind;
    double obs, *pot, *NLL, *g, *nodePot, *edgePot, *Xnode, *Xedge, Z;
    
    /* Input */
    NLL = mxGetPr(prhs[0]);
    g = mxGetPr(prhs[1]);
    i = (int)mxGetScalar(prhs[2])-1;
    nodePot = mxGetPr(prhs[3]);
    edgePot = mxGetPr(prhs[4]);
    edgeEnds = (int*)mxGetPr(prhs[5]);
    V = (int*)mxGetPr(prhs[6]);
    E = (int*)mxGetPr(prhs[7]);
    nStates = (int*)mxGetPr(prhs[8]);
    nodeMap = (int*)mxGetPr(prhs[9]);
    edgeMap = (int*)mxGetPr(prhs[10]);
    Xnode = mxGetPr(prhs[11]);
    Xedge = mxGetPr(prhs[12]);
    Y = (int*)mxGetPr(prhs[13]);
	
	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32")||!mxIsClass(prhs[7],"int32")||!mxIsClass(prhs[8],"int32")||!mxIsClass(prhs[9],"int32")||!mxIsClass(prhs[10],"int32")||!mxIsClass(prhs[13],"int32"))
		mexErrMsgTxt("i, edgeEnds, V, E, nStates, nodeMap, edgeMap, Y must be int32");
        
    /* Compute Sizes */
    nInstances = mxGetDimensions(prhs[11])[0];
    nNodeFeatures = mxGetDimensions(prhs[11])[1];
    nNodes = mxGetDimensions(prhs[3])[0];
    nEdgeFeatures = mxGetDimensions(prhs[12])[1];
    nEdges = mxGetDimensions(prhs[5])[0];
    maxState = getMaxState(nStates, nNodes);
	
    /* Allocate */
    pot = mxCalloc(maxState,sizeof(double));
    
    for(n = 0; n < nNodes; n++) {
        
        /* Compute Potential of each state with neighbors fixed */
        for(s = 0; s < nStates[n]; s++) {
            pot[s] = nodePot[n + nNodes*s];
        }
		
        for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++) {
            e = E[Vind]-1;
            n1 = edgeEnds[e]-1;
            n2 = edgeEnds[e+nEdges]-1;
            
            for(s = 0; s < nStates[n]; s++) {
                if(n == n1) {
                    pot[s] *= edgePot[s + maxState*(Y[i + nInstances*n2]-1 + maxState*e)];
                }
                else {
                    pot[s] *= edgePot[Y[i + nInstances*n1]-1 + maxState*(s + maxState*e)];
                }
            }
            
        }
        
        /* Local Normalizing Constant */
        Z = 0;
        for(s = 0; s < nStates[n]; s++) {
            Z += pot[s];
        }
        
        /* Update Objective */
        *NLL -= log(pot[Y[i + nInstances*n]-1]);
        *NLL += log(Z);
     
        /* Normalize */
        for(s = 0; s < nStates[n]; s++) {
            pot[s] /= Z;
        }
        
        /* Update gradient of node parameters */
        for(s = 0; s < nStates[n]; s++) {
            for(f = 0; f < nNodeFeatures; f++) {
                if(nodeMap[n + nNodes*(s + maxState*f)] > 0) {
                    if(s == Y[i + nInstances*n]-1)
                        obs = 1;
                    else 
                        obs = 0;
                    g[nodeMap[n + nNodes*(s + maxState*f)]-1] += Xnode[i + nInstances*(f + nNodeFeatures*n)]*(pot[s]-obs);
                }
            }
        }
        
        /* Update gradient of edge parameters */
        for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++) {
            e = E[Vind]-1;
            n1 = edgeEnds[e]-1;
            n2 = edgeEnds[e+nEdges]-1;
         
            for(s = 0; s < nStates[n]; s++) {
                if(n == n1) {
                    s1 = s;
                    s2 = Y[i + nInstances*n2]-1;
                }
                else {
                    s2 = s;
                    s1 = Y[i + nInstances*n1]-1;
                }
                for(f = 0; f < nEdgeFeatures; f++) {
                    if(edgeMap[s1 + maxState*(s2 + maxState*(e + nEdges*f))] > 0) {
                        if(s == Y[i + nInstances*n]-1)
                            obs = 1;
                        else
                            obs = 0;
                        g[edgeMap[s1 + maxState*(s2 + maxState*(e + nEdges*f))]-1] += Xedge[i + nInstances*(f + nEdgeFeatures*e)]*(pot[s]-obs);
                    }
                }
            }
        }
    }
    
    mxFree(pot);
        
}
