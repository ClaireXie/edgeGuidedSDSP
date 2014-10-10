EdgeGuided SDSP
=======================
Edge Guided Single Depth Image Super-resolution

This code implements the approach for super resolution with a single depth image input in the follwoing paper:

http://staff.washington.edu/junx/publication/icip2014_edgeGuided.pdf

If you intend to use the source code, please cite the paper as follows:
J. Xie, R. Feris and M.T. Sun, "Edge Guided Single Depth Image Super Resolution," ICIP 2014


How to Use the Code
=======================
Dependences:
The code depends on UGM toolbox (http://www.cs.ubc.ca/~schmidtm/Software/UGM.html), which is included in the utils/ folder. Please also cite its related work in order to use the toolbox. Remember to compile the toolbox before going to the next step. 

1. Run compileFiles.m to compile all the necessary mex files.

2. Then simply run demoFramework.m. An examplar learned dictionary file is included in dictionaries/. We use the training images from http://visual.cs.ucl.ac.uk/pubs/depthSuperRes/. 

3. Some example input depth images are included in inputs/

4. If you want to do the training on your own, with the collected image data, use edgeScript.m


Depth Super-resolution Results
=======================
Please note that the quantitative evaluation result in this version is slightly lower than that was reported in the paper since we utilize some code optimization for efficiency concerns in this version. 

Please feel free to contact junx@uw.edu for any questions or bug reports.

