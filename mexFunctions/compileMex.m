clear; clc;

% compile the mex files
mex -largeArrayDims bilateralUpsample.cpp fibheap.cpp
mex -largeArrayDims createGraph.cpp
