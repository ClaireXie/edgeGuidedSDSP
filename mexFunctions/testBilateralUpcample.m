% test the mex functions

clc; clear; close all;
compileMex;

addpath('../mainCode/');
addpath('../funcs/');

load('tempResult.mat');



depthMedium = imresize(input, scale, 'bicubic');
window = 7;
sigma_d = 0.5;

window_half = floor(window/2);
template = computeDistTemp(scale*2*window_half+1);

[structA, structC0] = constructGraphEdge(scale*2*window_half+1);
structA = structA+transpose(structA);
structC = structC0+transpose(structC0);
G =sparse(structC);
 
NHOOD = [0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);
edgeMap = edges+imdilate(edges, se)+1;
%w1 = constructGraph(G, double(edgeMap), nnz(G), n);


% change from matlab index to c++ index (except template)
edges = edges';
depthMedium = depthMedium';
edgeMap = edgeMap';
input = input';
offset = -1;

[highRes] = bilateralUpsample(input, depthMedium, double(edges), template, ... 
    window, sigma_d, scale, offset, G, nnz(G), edgeMap);

% convert back to matlab index
highRes = (reshape(highRes, size(edges, 1), size(edges, 2)))';
figure; imshow(uint8(highRes));

