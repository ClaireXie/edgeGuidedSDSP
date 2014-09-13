#include <math.h> /* Needed for the ceil() prototype */
#include "mex.h"
#include <cstdio>
#include <cmath>

#include <iostream>
#include <fstream>
#include <vector>

#include "fibheap.h"

using namespace std;



#define PI 3.1415926

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

void extractShortestPath(double *temp, vector< vector<Point> > &list, int m, int n) {

	// c++ index for k
	for (int k = 0; k < m*n; k++) {

		// matlab index for the followings
		vector<Point> path;
		for (int j = 0; j < m; j++) {
			for (int i = 0; i < n; i++) {
				int index =m*n*k + j*n + i;
				if (temp[index]!=0) {
					path.push_back(Point(j, i));
				}
			}
		}
		list.push_back(path);
	}
}

float normpdf(float value, float mean, float sigma) {
	return 1.0/(sigma*sqrt(2*PI))*exp(-(value-mean)*(value-mean)/(2.0*sigma*sigma));
}


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

// convert predecessor indices to shortest paths
void pred2path(vector <int> P, int s, int t, vector<int> &path) {
    // c++ index
    int ti = t;
    while (ti != -1 && ti != s) {
        path.push_back(ti);
        ti = P[ti];
    }
    path.push_back(s);
}

void dilatePath(vector <int> &path, int n) {
	
	// put path into an indexing map
	bool map[n*n]; 
	int m = path.size();

	for (int i = 0; i<n*n; i++) {
		map[i] = false;
	}

	for (int i = 0; i < m; i++) {
		map[path[i]] = true;
	}

	// do not dilate the starting and ending point
	for (int i = 1; i < m-1; i++) {
		int y = path[i]/n;
	    int x = path[i]%n;

	    // x+1, x-1, y+1 y-1
	    // don't consider ordering
	    int index;
	    if (y > 0) {
	    	index = n*(y-1) + x;
	    	if (!map[index]) { path.push_back(index); map[index] = true;}
	    }
	    if (x > 0) {
	    	index = n*(y) + x-1;
	    	if (!map[index]) {path.push_back(index); map[index] = true;}
	    }
	    if (y < n-1) {
	    	index = n*(y+1) + x;
	    	if (!map[index]) {path.push_back(index); map[index] = true;}
	    }
	    if (x < n-1) {
	    	index = n*(y) + x+1;
	    	if (!map[index]) {path.push_back(index); map[index] = true;}
	    }

	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
	// input: 
	// low depth + medum depth (formed by bicubic interpolation)
	// guided edges (dilted)
	// shortest path template
	// window_size
	// sigma_d
	// m and n: row and column
	// scale (int)
	// offset

	if (nrhs != 11) {
		mexErrMsgTxt("11 inputs required: \n low_depth medium_depth(using bicubic interp) edges shortestPath_template window_size sigma_d scale offset \n");
	}

	if (nlhs != 1) {
		mexErrMsgTxt("1 output required: \n highRes_depth \n");
	}

	//load images (be careful about c++ index and matlab index)
	double *lowDepth = (double *) mxGetData(prhs[0]);
	double *mediumDepth = (double *) mxGetData(prhs[1]);
	double *edge = (double *) mxGetData(prhs[2]);
	int n = mxGetN(prhs[2]);     // rows
	int m = mxGetM(prhs[2]);     // cols
	double *temp = (double *) mxGetData(prhs[3]);

	// load scalars
	float window_size = mxGetScalar(prhs[4]);
	float sigma_d = mxGetScalar(prhs[5]);

	int scale = mxGetScalar(prhs[6]);
	int offset = mxGetScalar(prhs[7]);
	int nzero = mxGetScalar(prhs[9]);

	double *edgeDilated = (double *) mxGetData(prhs[10]);


	plhs[0] = mxCreateNumericMatrix(n, m, mxDOUBLE_CLASS, mxREAL);
	double *depthHigh = (double *)mxGetData(plhs[0]);


	int window_half = window_size/2;
	int l = (scale*2*window_half+1)*(scale*2*window_half+1);
	//plhs[1] = mxCreateSparse(l, l, nzero, mxREAL);
	mxArray *edgeWeight = mxCreateSparse(l, l, nzero, mxREAL);



	cout << "============================================"<<endl;
	cout << "Parameters: "<<endl;
	cout << "window_size = "<<window_size<<endl;
	cout << "sigma_d = "<<sigma_d<<endl;
	cout << "(m , n) = ("<<m<<", "<<n<<")"<<endl;
	cout << "scale = "<<scale<<endl;
	cout << "offset = "<<offset<<endl;
	cout << "nzero = "<<nzero<<endl;

	int i, j;

	vector< vector<Point> > list;
	extractShortestPath(temp, list, 2*window_half*scale+1, 2*window_half*scale+1);
	
	// use c++ index 
	for (int idx = 0; idx < m*n; idx++) {

		i = idx/m;
        j = idx%m;

        if (i < scale*window_half || i >= n-scale*window_half-offset ||
            j < scale*window_half || j >= m-scale*window_half-offset) 
            continue;

        float avg = 0;
        float normalize = 0;
        
        // determine if the patch contains edges
        int edgeValue = 0;
        bool containEdge = false;
        vector<int> P;
        vector<double> D;
        for (int iy = -scale*window_half; iy <=  scale*window_half; iy++ ) {
        	for (int ix = -scale*window_half; ix <=  scale*window_half; ix++ ) {

        		int x = j + ix;
        		int y = i + iy;
        		int index = x + y*m;
        		edgeValue += edge[index];
        	}
        }

        // if this patch does not contain edge
        
        if (edgeValue > 0) 
        {
        	containEdge = true;
        	 // compute graph below
        	computeGraph(scale*2*window_half+1, edgeDilated, Point(j, i), prhs[8], edgeWeight, nzero, n);
        	// 0 base index
        	dijkstra( l, (l-1)/2, D, P, edgeWeight);
        }

        for (int iy = -scale*window_half; iy <=  scale*window_half; iy++ ) {
        	for (int ix = -scale*window_half; ix <=  scale*window_half; ix++ ) {
        		int x = j + ix;
        		int y = i + iy;

        		if ((x+offset)%scale == 0 && (y+offset)%scale == 0) {
        			int index = (iy+scale*window_half)*(scale*2*window_half+1) + (ix+scale*window_half);

        			bool crossEdge = false;
        			int crossEdgeCount = 0;

        			if (containEdge) {

        				// exam the shortest path from the template
	        			vector <int> path;
	        			pred2path(P, (l-1)/2, index, path);
	        			dilatePath(path, scale*2*window_half+1);

	        			// traverse the shortest path
	        			// also dilate the path

	        			for (int k = 0; k<path.size(); k++) {
	        				int yy = path[k]/(scale*window_half*2+1) + i - scale*window_half;
	        				int xx = path[k]%(scale*window_half*2+1) + j - scale*window_half;

	        				//cout << "A " << xx <<" "<<yy<<endl;

	        				if (edge[xx+ yy*m] != 0) {
	        					crossEdge = true;
	        					crossEdgeCount ++;;
	        				}
	        			}
	        		}

	        		// if (path.size() != list[index].size()) { 
	        			// 	cout << int(path.size()) - int(list[index].size())<<endl;
	        			// }

	        		/*else {

	        			// using template
	        			for (int k = 0; k<list[index].size(); k++) {
	        				int xx = list[index][k].x + j - scale*window_half;
	        				int yy = list[index][k].y + i - scale*window_half;

	        				//cout << "B " << xx <<" "<<yy<<endl;

	        				// edge should be in C++ index
	        				if (edge[xx+ yy*m] != 0) {
	        					crossEdge = true;
	        					crossEdgeCount ++;;
	        				}
	        			}
        			}*/

        			if (!crossEdge || (edge[j+ i*m] == 1 && crossEdgeCount == 1))
        			{
        				// perform bilateral filtering
        				float value = lowDepth[(x+offset)/scale + (y+offset)/scale * m/scale];
        				float d = sqrt((float(i)/scale-float(y)/scale)*(float(i)/scale-float(y)/scale)
        					+(float(j)/scale-float(x)/scale)*(float(j)/scale-float(x)/scale)); 
        				float g_d = normpdf(d, 0, sigma_d);
        				avg += value*g_d;
        				normalize += g_d;
        			}


        		}
        	}
        }

        if (normalize > 0) {
        	depthHigh[idx] = avg/normalize;
        }
        else {
        	// need to change from matlab index to c++ index
        	depthHigh[idx] = mediumDepth[idx];
        }
	}
}

