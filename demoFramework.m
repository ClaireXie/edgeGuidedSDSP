
% Demo of the edge guided single depth image super resolution + evaluation
% (c)2016 Jun Xie

clc; clear; close all;

addpath('mainCode/');
addpath('funcs/');
addpath('funcs/ann/');

% list of the input testing images
names{1} = 'cones';
names{2} = 'teddy';
names{3} = 'tsukuba';
names{4} = 'venus';
names{5} = '11_250';
names{6} = 'scan021_cave_statue_smNN';
names{7} = 'scan030_cave_statue_smNN';
names{8} = 'scan042_cave_statue_smNN';

% parameters to change
%----------------%
scaleFact = 4;
testIndex = 1;
inputFile = names{testIndex};
%----------------%
global window sigma_d;

% parameters setting
%----------------------------------%
w1 = 3;      % unary weights -- weights for the shock edge map
w2 = 1;      % pairwise weights
localSize = 1;  
window = 7;
sigma_d = 0.5;
scale = 4;
threshold = 0.08;
show = 1;
self_similarity = 1;

fprintf(['runnning image ',inputFile,'\n']);
disp('=======================================');

[highres,edges] = mrfLearning(names, testIndex, w1, w2, localSize, scale, ... 
    threshold, self_similarity, show);

if (testIndex ~= 5)
    % run evaluation
    border = 2*window;
    runEvaluation(inputFile, scale, highres, edges, scaleFact, border, threshold, ... 
        0, testIndex);
end
