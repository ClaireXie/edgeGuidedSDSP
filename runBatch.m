% the total script including:
% 1) super resolution
% 2) run evalution

clc;clear;close all;


if matlabpool('size') == 0
    matlabpool open local 4
end


addpath('mainCode/');

names=cell(4);
names{1}='cones';
names{2}='teddy';
names{3}='tsukuba';
names{4}='venus';
names{5}='11_250';
names{6}='scan042_cave_statue_smNN';
names{7}='scan030_cave_statue_smNN';
names{8}='scan021_cave_statue_smNN';

global window sigma_d sigma_c;
window = 7;
sigma_d = 0.5;
sigma_c = 0.05;   %10 or 0.05
w1 = 3;
w2 = 1;
localSize = 1;  %1
show = 1;
scale = 4;
threshold = [0.08, 0.1, 0.1, 0.1 0.05, 0.06, 0.06, 0.06];

para.window = window;
para.sigma_d = sigma_d;
para.sigma_c = sigma_c;
para.w1 = w1;
para.w2 = w2;
para.localSize = localSize;
para.scale = scale;


% run the code for each image
for i = 1:4
    
    if (i == 1 || i == 2)
        scaleFact = 4;
    elseif (i == 3)
        scaleFact = 16;
    elseif (i == 4)
        scaleFact = 8;
    end
    
    para.threshold = threshold(i);
    
    inputFile = names{i};
    fprintf(['runnning image ',names{i},'\n']);
    [highres{i}, edges{i}] = mrfLearning2(names, i, w1, w2, localSize, scale, threshold(i), show);
    
    border = 2*window;
    offset = 0;
    runEvaluation(inputFile, scale, highres{i}, edges{i}, scaleFact, border, threshold(i), 0);
    
    highRes = highres{i};
    highEdges = edges{i};
    save(sprintf('outputs/newExp/%s', inputFile), 'highRes', 'highEdges', 'para');
end
