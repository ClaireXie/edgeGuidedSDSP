EdgeGuided DSP
=======================
Edge Guided Single Depth Image Super-resolution

This code implements the approach for super resolution with a single depth image input in the follwoing paper:

http://staff.washington.edu/junx/publication/icip2014_edgeGuided.pdf

J. Xie, R. Feris and M.T. Sun, "Edge Guided Single Depth Image Super Resolution," ICIP 2014

How to Use the Code
=======================
1. Run compileFiles.m to compile all the necessary mex files.

2. Then simply run demoFramework.m. An examplar learned dictionary file is included in dictionaries/. We use the training images from http://visual.cs.ucl.ac.uk/pubs/depthSuperRes/.

3. If you wanna do the training on your own, with the collected image data, run edgeScript.m

MISC
=======================
PLease note that the result in this version of code is slightly lower than that was reported in the paper since we utilize some code optimization for efficiency concerns in this version. 

PLease feel free to contact junx@uw.edu for if you have questions regrading to the code. 

