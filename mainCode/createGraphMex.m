function [ unary edgePots edgeEnds] = createGraphMex(candidateH, candidateHTrans, index, diff, edgeMap)

addpath('mexFunctions/');

fprintf('construct the edge weights...\n');

% define the unary term
unary = diff.^2;

patchSize = sqrt(size(candidateH,3));
[n, m] = size(edgeMap);
labels = size(candidateH,2);

candidateHTrans = permute(candidateHTrans, [1,3,2]);
[structure, edgePots, edgeEnds] = createGraph(index, candidateHTrans, ... 
    patchSize, m, n, labels, size(index, 1));
edgePots = reshape(edgePots, [size(edgePots,1), labels, labels]);