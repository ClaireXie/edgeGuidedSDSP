#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
    int n, s,done,maxInd,e,n1,n2,Vind,s1,s2,
    nNodes, nEdges, maxState, dims[3],
    *edgeEnds, *nStates, *V, *E,*y,*yMax;
    
    double *pot,maxVal,
    *nodePot, *edgePot;
    
   /* Input */
    
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    edgeEnds = (int*)mxGetPr(prhs[2]);
    nStates = (int*)mxGetPr(prhs[3]);
    V = (int*)mxGetPr(prhs[4]);
    E = (int*)mxGetPr(prhs[5]);
    y = (int*)mxGetPr(prhs[6]);
	
	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32"))
        mexErrMsgTxt("edgeEnds, nStates, V, E, y must be int32");
    
   /* Compute Sizes */
    
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[2])[0];

   /* Output */
    dims[0] = nNodes;
    dims[1] = 1;
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    yMax = mxGetData(plhs[0]);
    
   /* Initialize */
    pot = mxCalloc(maxState,sizeof(double));
    
   /* Start at Maximum NodePot */
    /*for(n=0;n<nNodes;n++)
    {
        maxVal = -1;
        maxVal = -1;
        for(s=0;s<nStates[n];s++)
        {
            if(nodePot[n + nNodes*s] > maxVal)
            {
                maxVal = nodePot[n + nNodes*s];
                maxInd = s;
            }
        }
        y[n] = maxInd;
    }*/
    
	for(n = 0; n < nNodes; n++)
    {
        yMax[n] = y[n]-1;
    }
	
    done = 0;
    while(!done)
    {
        done = 1;
        for(n = 0; n < nNodes; n++)
        {
           /* Compute Node Potential */
            for(s = 0; s < nStates[n]; s++)
            {
                pot[s] = nodePot[n + nNodes*s];
            }
            
           /* Iterate over Neighbors */
            for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
            {
                e = E[Vind]-1;
                n1 = edgeEnds[e]-1;
                n2 = edgeEnds[e+nEdges]-1;
                 
                /* Multiply Edge Potentials */
                if(n == n1)
                {
                   for(s = 0; s < nStates[n]; s++)
                   {
                        pot[s] *= edgePot[s+maxState*(yMax[n2] + maxState*e)];
                   }
                    
                }
                else
                {
                    for(s = 0; s < nStates[n]; s++)
                    {
                        pot[s] *= edgePot[yMax[n1]+maxState*(s + maxState*e)];
                    }
                }
                
            }
            
            
            
           /* Assign to Maximum State */
            maxVal = -1;
            for(s=0;s<nStates[n];s++)
            {
                if(pot[s] > maxVal)
                {
                    maxVal = pot[s];
                    maxInd = s;
                }
            }
            if (maxInd != yMax[n])
            {
                yMax[n] = maxInd;
                done = 0;
            }
            
        }
    }
    
    for(n = 0; n < nNodes; n++)
    {
        yMax[n] = yMax[n]+1;
    }
    
    
   /* Free memory */
    mxFree(pot);
}
