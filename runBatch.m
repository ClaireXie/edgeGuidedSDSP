
% Demo of the edge guided single depth image super resolution + evaluation
% (batched version)
% (c)2016 Jun Xie

clc;clear;close all;
warning off;

addpath('mainCode/');

names{1}='cones';
names{2}='teddy';
names{3}='tsukuba';
names{4}='venus';
names{5}='11_250';
names{6}='scan042_cave_statue_smNN';
names{7}='scan030_cave_statue_smNN';
names{8}='scan021_cave_statue_smNN';

global window sigma_d;
window = 7;
sigma_d = 0.5;
w1 = 3;
w2 = 1;
localSize = 1;
scale = 4;
threshold = [0.08, 0.1, 0.1, 0.1 0.05, 0.06, 0.06, 0.06];

para.window = window;
para.sigma_d = sigma_d;
para.w1 = w1;
para.w2 = w2;
para.localSize = localSize;
para.scale = scale;

% flags:
print2File = 0;

% ---------- Enable self-similarity --------------%
self_similarity = 1;
% ------------------------------------------------%

% ---------- Enable visualize results ------------%
show = 0;
% ------------------------------------------------%

% run the code for each image
for i = 1:numel(names)
    
    % The scale fact is fixed (please refer to the Middlebury Dataset)
    scaleFact = 1;
    if (i == 1 || i == 2)
        scaleFact = 4;
    elseif (i == 3)
        scaleFact = 16;
    elseif (i == 4)
        scaleFact = 8;
    end
    
    para.threshold = threshold(i);
    
    inputFile = names{i};
    fprintf('runnning image %s X %d \n', names{i}, scale);
    [highres{i}, edges{i}] = mrfLearning(names, i, w1, w2, localSize, ... 
        scale, threshold(i), self_similarity, show);
    
    border = 2*window;
    if (i ~= 5)
        runEvaluation(inputFile, scale, highres{i}, scaleFact, border, print2File, i);
    end
    
    highRes = highres{i};
    highEdges = edges{i};
    save(sprintf('outputs/%s', inputFile), 'highRes', 'highEdges', 'para');
end
