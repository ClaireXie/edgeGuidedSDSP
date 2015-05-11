EdgeGuided SDSP
=======================
Edge Guided Single Depth Image Super-resolution

This code implements the approach for super resolution with a single depth image input in the following paper:

http://public.genius-field.com/junx/publication/icip2014_edgeGuided.pdf

If you intend to use the source code, please cite the paper as follows:
J. Xie, R. Feris and M.T. Sun, "Edge Guided Single Depth Image Super Resolution," ICIP 2014


Dependences
=======================
1. UGM toolbox (http://www.cs.ubc.ca/~schmidtm/Software/UGM.html), included in the utils/ folder. Please also cite its related work in order to use the toolbox. Remember to compile the toolbox before going to the next step. 


2. ANN (https://www.cs.umd.edu/~mount/ANN/), *optional* for a more efficient k-nearest neighbor search implementation. 
The source code of ANN is included in funs/ANN/. We have provided with a Matlab wrapper for ANN library. 
This feature is enabled by default. To disable it, simply set line68@mainCode/mrfLearning.m: 
	useANN = 1;
to
	useANN = 0;

Then run make.m to compile the mex file of the ANN wrapper. 


How to Use the Code
=======================
1. Run compileFiles.m to compile all the necessary mex files.

2. Then simply run demoFramework.m. An examplar learned dictionary file is included in dictionaries/. We use the training images from http://visual.cs.ucl.ac.uk/pubs/depthSuperRes/. 

3. Some example input depth images are included in inputs/

4. If you want to do the training on your own, with the collected image data, use edgeScript.m

Note: The code has been tested under 64bit Linux and Windows platform with Matlab 2013b and 2014b installed. 


Depth Super-resolution Results
=======================
Please note that the quantitative evaluation result in this version is slightly lower than that was reported in the paper since we utilize some code optimization for efficiency concerns in this version. 

Please feel free to contact junx@uw.edu for any questions or bug reports.

