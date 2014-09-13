clc; clear; close all;


addpath('mainCode/');

% parameters to change
%----------------%
names{1} ='cones';
inputFile = names{1};
scaleFact = 4;
%----------------%
global window sigma_d;

% multi-thread
if matlabpool('size') == 0
    matlabpool open local 4
end

% parameters setting
%----------------------------------%
w1 = 3;      %unary weights -- weights for the shock edge map
              %              -- weights for the original edge map = 1
w2 = 1;       %pairwise weights
localSize = 1;  
window = 7;
sigma_d = 0.5;
scale = 4;
threshold = 0.08;

show = 1;


fprintf(['runnning image ',inputFile,'\n']);
disp('=======================================');

totalTime = tic;
[highres,edges] = mrfLearning2(names, 1, w1, w2, localSize, scale, threshold, show);
fprintf('Total time = \n');
toc(totalTime);

% with only knn (unary term)
%[highres,edges] = mrfLearning(names, 1, w1, w2, localSize, scale, threshold);

% run evaluation

border = 2*window;
runEvaluation(inputFile, scale, highres, edges, scaleFact, border, threshold,0);


