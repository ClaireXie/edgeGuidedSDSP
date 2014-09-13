#include <math.h>
#include "mex.h"
#include "UGM_common.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
    int n, s,e,n1,n2,neigh,Vind,s1,s2,
    nNodes, nEdges, maxState, dims[3],
    iter,maxIter,
    *edgeEnds, *nStates, *V, *E,*y;
    
    double *nodePot, *edgePot, *nodeBel, *edgeBel, *logZ,
    *oldNodeBel,*b,z,U1,U2,S1;
    
   /* Input */
    
    nodePot = mxGetPr(prhs[0]);
    edgePot = mxGetPr(prhs[1]);
    edgeEnds = (int*)mxGetPr(prhs[2]);
    nStates = (int*)mxGetPr(prhs[3]);
    V = (int*)mxGetPr(prhs[4]);
    E = (int*)mxGetPr(prhs[5]);
    maxIter = ((int*)mxGetPr(prhs[6]))[0];
    
	if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32")||!mxIsClass(prhs[4],"int32")||!mxIsClass(prhs[5],"int32")||!mxIsClass(prhs[6],"int32"))
		mexErrMsgTxt("edgeEnds, nStates, V, E, maxIter must be int32");
	
   /* Compute Sizes */
    
    nNodes = mxGetDimensions(prhs[0])[0];
    maxState = mxGetDimensions(prhs[0])[1];
    nEdges = mxGetDimensions(prhs[2])[0];
    
   /* Output */
    plhs[0] = mxCreateDoubleMatrix(nNodes,maxState,mxREAL);
    dims[0] = maxState;
    dims[1] = maxState;
    dims[2] = nEdges;
    plhs[1] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL);
    nodeBel = mxGetPr(plhs[0]);
    edgeBel = mxGetPr(plhs[1]);
    logZ = mxGetPr(plhs[2]);
    
    /* Initialize */
    for(n = 0; n < nNodes; n++)
    {
        z = 0;
        /* Initialize nodeBel to nodePot */
        for(s = 0;s < nStates[n];s++)
        {
            nodeBel[n + nNodes*s] = nodePot[n + nNodes*s];
            z += nodeBel[n + nNodes*s];
        }
         /* Normalize nodeBel */
        for(s = 0; s < nStates[n];s++)
        {
            nodeBel[n + nNodes*s] /= z;
        }
    }
    
    oldNodeBel = mxCalloc(nNodes*maxState,sizeof(double));
    b = mxCalloc(maxState,sizeof(double));
    
    for(iter = 0; iter < maxIter; iter++)
    {
        
      /* oldNodeBel = nodeBel */
        for(n=0;n<nNodes;n++)
        {
            for(s=0;s<nStates[n];s++)
            {
                oldNodeBel[n+nNodes*s] = nodeBel[n+nNodes*s];
            }
        }
        
        for(n=0;n<nNodes;n++)
        {
            
            for(s = 0; s < nStates[n]; s++)
                b[s] = 0;
            
            /* Collect Messages from all Neighbors */
            for(Vind = V[n]-1; Vind < V[n+1]-1; Vind++)
            {
                e = E[Vind]-1;
                n1 = edgeEnds[e]-1;
                n2 = edgeEnds[e+nEdges]-1;
                
                if(n == n2)
                {
                    neigh = n1;
                    for(s = 0; s < nStates[n]; s++)
                    {
                        for(s2 = 0; s2 < nStates[neigh]; s2++)
                        {
                            b[s] += nodeBel[neigh + nNodes*s2]*log(edgePot[s2 + maxState*(s + maxState*e)]);
                        }
                    }
                }
                else
                {
                    neigh = n2;
                    for(s = 0; s < nStates[n]; s++)
                    {
                        for(s2 = 0; s2 < nStates[neigh]; s2++)
                        {
                            b[s] += nodeBel[neigh + nNodes*s2]*log(edgePot[s + maxState*(s2 + maxState*e)]);
                        }
                    }
                }
                
            }
            
            /* Update Local Belief */
            z = 0;
            for(s = 0; s < nStates[n]; s++)
            {
                b[s] = nodePot[n + nNodes*s]*exp(b[s]);
                z += b[s];
            }
            for(s = 0; s < nStates[n]; s++)
            {
                nodeBel[n + nNodes*s] = b[s]/z;
            }
            
        }
        
        
        /* if sum(abs(nodeBel(:)-oldNodeBel(:))) < 1e-4; break; */
        z = 0;
        for(n=0;n<nNodes;n++)
            for(s=0;s<nStates[n];s++)
                z += absDif(nodeBel[n+nNodes*s],oldNodeBel[n+nNodes*s]);
        
        if(z < 1e-4)
        {
            break;
        }
        
    }
    
    /*if(iter == maxIter)
    {
        printf("MF reached maxIter of %d iterations\n",maxIter);
    }*/
    /* printf("Stopped after %d iterations\n",i); */
    
    
    /* Compute edgeBel */
    for(e = 0; e < nEdges; e++)
    {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        for(s1 = 0; s1 < nStates[n1]; s1++)
        {
            for(s2 = 0; s2 < nStates[n2]; s2++)
            {
                edgeBel[s1+maxState*(s2+maxState*e)] = nodeBel[n1+nNodes*s1]*nodeBel[n2+nNodes*s2];
            }
        }
    }
    
    /* Compute Gibbs Free Energy */
    U1 = 0;
    U2 = 0;
    S1 = 0;
    
    for(n = 0; n < nNodes; n++)
    {
        for(s = 0; s < nStates[n]; s++)
        {
            U1 += nodeBel[n+nNodes*s]*log(nodePot[n+nNodes*s]);
            if(nodeBel[n + nNodes*s] > 1e-10)
            {
                S1 += nodeBel[n+nNodes*s]*log(nodeBel[n+nNodes*s]);
            }
        }
    }
    
    for(e = 0; e < nEdges; e++)
    {
        n1 = edgeEnds[e]-1;
        n2 = edgeEnds[e+nEdges]-1;
        
        for(s1 = 0; s1 < nStates[n1]; s1++)
        {
            for(s2 = 0; s2 < nStates[n2]; s2++)
            {
                U2 += nodeBel[n1+nNodes*s1]*nodeBel[n2+nNodes*s2]*log(edgePot[s1+maxState*(s2+maxState*e)]);
            }
        }
    }
    logZ[0] = U2 + U1 - S1;
    
   /* Free memory */
    mxFree(oldNodeBel);
    mxFree(b);
}
