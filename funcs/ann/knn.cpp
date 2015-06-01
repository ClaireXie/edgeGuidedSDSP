#include "ANN/ANN.h"
#include <vector>
#include <fstream>
#include <iostream>
#include <ctime>

using namespace std;

/* read data from file
 * @params ANNpointArray data structure (num*dim)
 * @params dimension of the data
 * @params name of the input file
 */ 
int readPt(ANNpointArray &dataPts, int dim, string filename)			
{
	ifstream infile;
	infile.open(filename.c_str());
	int count = 0;
	while (!infile.eof()) {
		ANNpoint p = annAllocPt(dim); 
		for (int i = 0; i<dim; i++) {
			double x;
			infile >> x;
			if (!infile.eof())
				p[i] = x;
		}
		if (!infile.eof()) {
			dataPts[count] = p;
			count ++;
		}
	}
	infile.close();

	return count;
}


int main (int argc, char **argv)
{
	ANNpointArray dataPts;
	ANNpointArray queryPt;
	vector<ANNidxArray>	nnIdx;	
	vector<ANNdistArray> dists;	
	ANNkd_tree*	kdTree;

	int dim = atoi(argv[3]);
	int k = atoi(argv[4]);

	int maxPts = 100000;
	dataPts = annAllocPts(maxPts, dim);	
	queryPt = annAllocPts(maxPts, dim);	

	int num = readPt(dataPts, dim, argv[1]);
	cout << "Load " << num << " Points. " << endl;

	int num_query = readPt(queryPt, dim, argv[2]);
	cout << "Load " << num_query << " Query points. " << endl;
	
	clock_t start = clock();

	// build the kd tree
	cout << "Building the kd tree..." << endl;
	kdTree = new ANNkd_tree(dataPts, num, dim);		

	// query
	cout << "knn searching..." << endl;
	for (int i = 0; i<num_query; i++) {
		ANNidxArray	nnIdx0 = new ANNidx[k];	
		ANNdistArray dists0 = new ANNdist[k];
		kdTree->annkSearch(	queryPt[i], k, nnIdx0, dists0, 0);				

		nnIdx.push_back(nnIdx0);
		dists.push_back(dists0);
	}

	cout << "Duration = " << ( std::clock() - start ) / (double) CLOCKS_PER_SEC << " sec." << endl;

	//print out index
	// for (int j = 0; j<num_query; j++) {
	// 	cout << j+1 << " ";
	// 	for (int i = 0; i < k; i++) {
	// 		cout << nnIdx[j][i]+1 << " ";
	// 	}
	// 	cout << endl;
	// }

	// free memory
	for (int i = 0; i<num_query; i++) {
		delete [] nnIdx[i];
    	delete [] dists[i];
	} 

	annDeallocPts(dataPts);
    annDeallocPts(queryPt);

	 delete kdTree;
    annClose();

	return 1;
}