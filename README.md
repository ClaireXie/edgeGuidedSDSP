EdgeGuided SDSP
=======================
Edge Guided Single Depth Image Super-resolution

This code implements the approach for super resolution with a single depth image input in this [paper](http://www.clairexie.org/resources/TIP16.pdf)

If you intend to use the source code, please cite the paper as follows:

```
Jun Xie, Rogerio. S. Feris and Ming-Ting Sun, "Edge-Guided Single Depth Image Super Resolution," 
in IEEE Transactions on Image Processing, vol. 25, no. 1, pp. 428-438, Jan. 2016.
```

Dependences
=======================
1. UGM toolbox (http://www.cs.ubc.ca/~schmidtm/Software/UGM.html), included in the utils/ folder. Please also cite its related work in order to use the toolbox. Remember to compile the toolbox before going to the next step. 


2. ANN (https://www.cs.umd.edu/~mount/ANN/), *optional* for a more efficient k-nearest neighbor search implementation. 
The source code of ANN is included in funs/ANN/. We have provided with a Matlab wrapper for ANN library. 
This feature is enabled by default. To disable it, simply set the flag in **mainCode/mrfLearning.m** from
	
	```
	useANN = 1; => useANN = 0;
	```

3. Download the trained dictionary:

	(If you only intend to run the self-similarity part, just ignore this step.)

	[Dictionary with upscaling factor = 3](http://www.clairexie.org/data/dictionaries/patchData_3_high.mat) [~167MB]

	[Dictionary with upscaling factor = 4](http://www.clairexie.org/data/dictionaries/patchData_4_high.mat) [~75MB]


How to Use the Code
=======================
1. Run compileFiles.m to compile all the necessary mex files. (Pre-built mex for Windows and Linux are included)

2. *demoFramework.m* is a simple demo script. *runBatch.m* is the batch script to run a couple of images.

4. Some example input depth images are included in inputs/

5. If you intend to do the training on your own, with the collected image data, use trainingScript.m. 
   We use the training images from http://visual.cs.ucl.ac.uk/pubs/depthSuperRes/. After getting the images, please modify the data directory in trainingScript.m.

6. You can switch to the self-similarity mode (without the training data), change from 

	```
	self_similarity = 1; => self_similarity = 0;
	```

	in **runBatch.m**/**demoFramework.m**

Note: The code has been tested under 64bit Linux and Windows platform with Matlab 2014b/2015a installed. 


Depth Super-resolution Results
=======================
Please note that the quantitative evaluation result in this version is slightly lower than that was reported in the paper since we utilize some code optimization for efficiency concerns in this version. 

Please feel free to contact xjsjtu88@gmail.com for any questions or bug reports.

Change Logs
=======================
2016-05-19: Add the self-similarity code, remove dictionary from repo to server


