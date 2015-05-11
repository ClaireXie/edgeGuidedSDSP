#include <math.h> /* Needed for the ceil() prototype */
#include "mex.h"
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>

#include "fibheap.h"

using namespace std;

class HeapNode : public FibHeapNode {
  double   N;
  long int IndexV;
  
public:
  HeapNode() : FibHeapNode() { N = 0; };
  virtual void operator =(FibHeapNode& RHS);
  virtual int  operator ==(FibHeapNode& RHS);
  virtual int  operator <(FibHeapNode& RHS);
  virtual void operator =(double NewKeyVal );
  virtual void Print();
  double GetKeyValue() { return N; };
  void SetKeyValue(double n) { N = n; };
  long int GetIndexValue() { return IndexV; };
  void SetIndexValue( long int v) { IndexV = v; };
};

void HeapNode::Print() {
  FibHeapNode::Print();
  mexPrintf( "%f (%d)" , N , IndexV );
}

void HeapNode::operator =(double NewKeyVal) {
  HeapNode Tmp;
  Tmp.N = N = NewKeyVal;
  FHN_Assign(Tmp);
}

void HeapNode::operator =(FibHeapNode& RHS) {
  FHN_Assign(RHS);
  N = ((HeapNode&) RHS).N;
}

int  HeapNode::operator ==(FibHeapNode& RHS) {
  if (FHN_Cmp(RHS)) return 0;
  return N == ((HeapNode&) RHS).N ? 1 : 0;
}

int  HeapNode::operator <(FibHeapNode& RHS) {
  int X;
  if ((X=FHN_Cmp(RHS)) != 0)
    return X < 0 ? 1 : 0;
    return N < ((HeapNode&) RHS).N ? 1 : 0;
};


class Point
{
    public:
    int x; 
    int y;
    Point(int x0, int y0) {x = x0; y = y0;}
    ~Point(){};
};

void dijkstra1( long int n, long int s, double *D1, double *P1, double *Gpr, mwIndex *Gir, mwIndex *Gjc) {
  int      finished;
  long int i, startInd, endInd, whichNeigh, nDone, closest;
  double   closestD, arcLength, INF, SMALL, oldDist;
  HeapNode *A, *hnMin, hnTmp; FibHeap *heap;
  INF=mxGetInf(); SMALL=mxGetEps();
  
  // setup heap
  if ((heap = new FibHeap) == NULL || (A = new HeapNode[n+1]) == NULL )
    mexErrMsgTxt( "Memory allocation failed-- ABORTING.\n" );
  heap->ClearHeapOwnership();
  
  // initialize
  for (i=0; i<n; i++) {
    if (i!=s) A[ i ] = (double) INF; else A[ i ] = (double) SMALL;
    if (i!=s) D1[ i ] = (double) INF; else D1[ i ] = (double) SMALL;
    P1[ i ] = -1;
    heap->Insert( &A[i] );
    A[ i ].SetIndexValue( (long int) i );
  }
  
  // Insert 0 then extract it, which will balance heap
  heap->Insert(&hnTmp); heap->ExtractMin();
  
  // loop over nonreached nodes
  finished = nDone = 0;
  while ((finished==0) && (nDone < n)) {
    hnMin = (HeapNode *) heap->ExtractMin();
    closest  = hnMin->GetIndexValue();
    closestD = hnMin->GetKeyValue();
    if ((closest<0) || (closest>=n))
      mexErrMsgTxt( "Minimum Index out of bound..." );
    D1[ closest ] = closestD;
    if (closestD == INF) finished=1; else {
      // relax all nodes adjacent to closest
      nDone++;
      startInd = Gjc[ closest   ];
      endInd   = Gjc[ closest+1 ] - 1;
      if( startInd!=endInd+1 )
        for( i=startInd; i<=endInd; i++ ) {
        whichNeigh = Gir[ i ];
        arcLength = Gpr[ i ];
        oldDist = D1[ whichNeigh ];
        if ( oldDist > ( closestD + arcLength )) {
          D1[ whichNeigh ] = closestD + arcLength;
          // c++ index
          P1[ whichNeigh ] = closest;
          hnTmp = A[ whichNeigh ];
          hnTmp.SetKeyValue( closestD + arcLength );
          heap->DecreaseKey( &A[ whichNeigh ], hnTmp );
        }
        }
    }
  }
  
  // cleanup
  delete heap; delete[] A;
}


// dijkastra shortest path algorithm
// nSrc = 1
void dijkstra( long int n, int sourcesIndx, vector<double> &D, vector<int> &P, const mxArray *G ) {
  // dealing with sparse array
    double *Gpr = mxGetPr(G);
    mwIndex *Gir = mxGetIr(G);
    mwIndex *Gjc = mxGetJc(G);
    int nSrc = 1;

    // allocate memory for single source results (automatically recycled)
    double *D1 = (double *) mxCalloc( n , sizeof( double ));
    double *P1 = (double *) mxCalloc( n , sizeof( double ));

    // loop over sources
    long int s, i, j;  
    i = 0;
    // run the dijkstra1 code for single source (0 indexed)
    s = sourcesIndx;
    if (s<0 || s > n-1) mexErrMsgTxt( "Source node(s) out of bound" );
    dijkstra1( n, s, D1, P1, Gpr, Gir, Gjc );

    // store results
    /*for( j=0; j<n; j++ ) {
      *( D + j*nSrc + i ) = *( D1 + j );
      *( P + j*nSrc + i ) = *( P1 + j );
    }*/

    for( j=0; j<n; j++ ) {
      D.push_back(*( D1 + j ));
      P.push_back(*( P1 + j ));
    }

}


void computeGraph(int n, double *dilateEdge, const Point center, const mxArray *G, mxArray *C, int nzero, int m)
{
	int half = (n-1)/2;
	// put outside the function
	//G = mxCreateSparse(m,n,nzmax,cmplx);
	mwIndex *irs,*jcs;
	double *sr;
	sr  = mxGetPr(G);
    irs = mxGetIr(G);
    jcs = mxGetJc(G);

    mwIndex *irs_c,*jcs_c;
    double *sr_c;
    sr_c  = mxGetPr(C);
    irs_c = mxGetIr(C);
    jcs_c = mxGetJc(C);

	int col = 0;
	int k = 0;
	Point offset(center.x-half, center.y-half);

    jcs_c[0] = 0;
    for (int i = 0; i<nzero; i++) {
    	int ind1 = irs[i];
    	int ind2 = col;
        irs_c[i] = irs[i];
    	k ++;
    	if (k >= jcs[col+1]) {
            jcs_c[col+1] = jcs[col+1];
            col ++;
    	}
    	int ind1_x = ind1%n;
    	int ind1_y = ind1/n;

    	int ind2_x = ind2%n;
    	int ind2_y = ind2/n;

        // dilatedEdge: matlab index
    	sr_c[i] = sr[i]*max(dilateEdge[(offset.y+ind1_y) + (offset.x+ind1_x)*m], 
            dilateEdge[offset.y+ind2_y + (offset.x+ind2_x)*m]);
    }
}

void pred2path(vector <int> P, int s, int t, vector<int> &path) {
    // c++ index
    int ti = t;
    while (ti != -1 && ti != s) {
        path.push_back(ti);
        ti = P[ti];
    }
    path.push_back(s);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {


  // if (nrhs < 4) {
  //   mexErrMsgTxt("4 inputs required: data_train (N1*dim), data_query(N2*dim), k, (verbose)");
  // }

  mwIndex *ir, *jc;
  double *s;
  ir = mxGetIr(prhs[0]);      /* Row indexing      */
  jc = mxGetJc(prhs[0]);      /* Column count      */
  s  = mxGetPr(prhs[0]);      /* Non-zero elements */ 

  double *edge = (double *) mxGetData(prhs[1]);
  int nzero = mxGetScalar(prhs[2]);
  int n = mxGetScalar(prhs[3]);
  int m = mxGetM(prhs[1]);

  cout<<"parameters:"<<endl;
  cout<<"m = "<<m<<endl;
  cout<<"n = "<<n<<endl;
  cout<<"nzero = "<<nzero<<endl;

  plhs[0] = mxCreateSparse(n*n, n*n, nzero, mxREAL);

  computeGraph(n, edge, Point(221, 179), prhs[0], plhs[0], nzero, m);
  int l = n*n;
  vector <int> P, path;
  vector <double> D;
  dijkstra( l, (n*n-1)/2, D, P, plhs[0]);

  pred2path(P, (n*n-1)/2, 3, path);

  /*for (int i = 0; i<path.size(); i++) {
      cout << path[i]<<endl;
  }*/
}