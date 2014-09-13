clc; clear; close all;
compileMex;

addpath('../mainCode/');
addpath('../funcs/');

load('tempResult.mat');

n = 5;
[structA, structC0] = constructGraphEdge(n);
structA = structA+transpose(structA);
structC = structC0+transpose(structC0);
G =sparse(structC);
 
NHOOD = [0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);
edgeMap = edges+imdilate(edges, se)+1;
w1 = constructGraph(G, double(edgeMap), nnz(G), n);

% for testing
i = 180;
j = 222;
half = (n-1)/2;
edgeTmp = edgeMap(i-half:i+half, j-half:j+half);   
[row, col] = find(structA == 1);
w2 = constructGraphWeight(row, col, edgeTmp, n, structC);

figure;imshow(full(w1));
figure;imshow(full(w2));

diff = full(w1)-w2;
figure;imshow(abs(diff));

[costs,paths] = dijkstra(structA,structC,(n*n+1)/2);
a = paths{4};