#include <math.h>
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
   /* Variables */
   int n, s,s1,s2,n1,n2,e, yInd,
   nNodes, nEdges, maxState, sizeYmap[3],
   *y,
   *edgeEnds, *nStates,*yMAP;
   
   double pot,maxPot,
   *nodePot, *edgePot;
   
   /* Input */
   
   nodePot = mxGetPr(prhs[0]);
   edgePot = mxGetPr(prhs[1]);
   edgeEnds = (int*)mxGetPr(prhs[2]);
   nStates = (int*)mxGetPr(prhs[3]);
   
   if (!mxIsClass(prhs[2],"int32")||!mxIsClass(prhs[3],"int32"))
        mexErrMsgTxt("edgeEnds and nStates must be int32");
   
   /* Compute Sizes */
   
   nNodes = mxGetDimensions(prhs[0])[0];
   maxState = mxGetDimensions(prhs[0])[1];
   nEdges = mxGetDimensions(prhs[2])[0];
   
   /* Output */
   sizeYmap[0] = nNodes;
   sizeYmap[1] = 1;
   plhs[0] = mxCreateNumericArray(2,sizeYmap,mxINT32_CLASS,mxREAL);
   yMAP = mxGetData(plhs[0]);
   
   /* Initialize */
   y = mxCalloc(nNodes,sizeof(int));
   maxPot = -1;
   while(1)
   {
      pot = 1;
      
   /* Node Potentials */
      for(n = 0; n < nNodes; n++)
      {
         pot *= nodePot[n + nNodes*y[n]];
      }
      
   /* Edge Potentials */
      for(e = 0; e < nEdges; e++)
      {
         n1 = edgeEnds[e]-1;
         n2 = edgeEnds[e+nEdges]-1;
         pot *= edgePot[y[n1] + maxState*(y[n2] + maxState*e)];
      }
      
      
      /* Compare potential of current configuration against best */
      if (pot > maxPot)
      {
        maxPot = pot;
        for(n = 0; n < nNodes;n++)
        {
            yMAP[n] = y[n]+1;
        }
      }
      
      
   /* Go to next y */
      
      for(yInd = 0; yInd < nNodes; yInd++)
      {
         y[yInd] += 1;
         if(y[yInd] < nStates[yInd])
         {
          break;  
         }
         else
         {
            y[yInd] = 0;
         }
      }

   /* Stop when we are done all y combinations */
      if(yInd == nNodes)
      {
         break;
      }
   }
   
   
   /* Free memory */
   mxFree(y);
}
