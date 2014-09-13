clc; clear; close all;

compileMex;

patchSize = 21;
n = 372; 
m = 448;
labels = 5;
load('tmpGraph');

candidateHTrans = permute(candidateHTrans, [1,3,2]);
[s, edgePots, edgeEnds] = createGraph(index, candidateHTrans, patchSize, m, n, labels, size(structure, 1));