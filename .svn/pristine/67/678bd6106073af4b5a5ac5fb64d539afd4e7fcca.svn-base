#include <math.h>
#include "mex.h"
#include "UGM_common.h"


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
    int i,n,s,e,n1,n2,Vind,maxIter,burnIn,
    nNodes, nEdges, maxState, dims[2],
    *edgeEnds, *nStates, *V, *E,*y,*S;
    
    double *pot,z,U,u,
    *nodePot, *edgePot;
    
   /* Input */
    
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    edgeEnds = (int*)mxGetPr(prhs[2]);
    nStates = (int*)mxGetPr(prhs[3]);
    V = (int*)mxGetPr(prhs[4]);
    E = (int*)mxGetPr(prhs[5]);
    maxIter = ((int*)mxGetPr(prhs[6]))[0];
    burnIn = ((int*)mxGetPr(prhs[7]))[0];
    y = (int*)mxGetPr(prhs[8]);
		
	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32")||!mxIsClass(prhs[7],"int32")||!mxIsClass(prhs[8],"int32"))
        mexErrMsgTxt("edgeEnds, nStates, V, E, maxIter, burnIn, y must be int32");
    
   /* Compute Sizes */
    
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[2])[0];
    
   /* Output */
    pot = mxCalloc(maxState,sizeof(double));
    dims[0] = nNodes;
    dims[1] = maxIter;
    plhs[0] = mxCreateNumericArray(2,dims,mxINT32_CLASS,mxREAL);
    S = mxGetData(plhs[0]);
    
   /* Initialize to States with highest node potentials*/
   /* for(n = 0; n < nNodes; n++)
    {
        u = -1;
        U = 0;
        for(s = 0; s < nStates[n]; s++)
        {
            if(nodePot[n+nNodes*s] > u)
            {
                u = nodePot[n+nNodes*s];
                U = s;
            }
        }
        y[n] = U;
    }
    */
	
	for(n = 0; n < nNodes; n++)
		y[n]--;
    
    for(i = 0; i < burnIn+maxIter; i++)
    {

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
                        pot[s] *= edgePot[s+maxState*(y[n2] + maxState*e)];
                   }
                    
                }
                else
                {
                    for(s = 0; s < nStates[n]; s++)
                    {
                        pot[s] *= edgePot[y[n1]+maxState*(s + maxState*e)];
                    }
                }
                
            }
            
            /* Normalize */
            z = 0;
            for(s = 0; s < nStates[n]; s++)
                z = z + pot[s];
            for(s = 0; s < nStates[n]; s++)
                pot[s] /= z;
            
            /* Display */
            for(s = 0; s < nStates[n]; s++)
            
            /* Sample Discrete State */
            U = rand()/((double)RAND_MAX + 1);
            u = 0;
            for(s = 0; s < nStates[n]; s++)
            {
                u += pot[s];
                if(u > U)
                {
                    break;
                }
            }
            y[n] = s;
        }
        
        if(i >= burnIn)
        {        /* Record Sample */
            for(n = 0; n < nNodes; n++)
            {
                S[n + nNodes*(i-burnIn)] = y[n]+1;
            }
        }
    }
    
   /* Free memory */
    mxFree(pot);
}
