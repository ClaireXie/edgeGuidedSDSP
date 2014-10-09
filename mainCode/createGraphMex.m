
function [ unary, edgePots, edgeEnds] = createGraphMex(candidateH, ... 
    candidateHTrans, index, diff, edgeMap)
% Mexed function for consturcting an MRF graph
%   Input: 
%        candidateH: high-res candidate patches
%        candidateH: high-res candidate patches
%        candidateHTrans: high-res candidate patches(distance transformed)
%        index: lookup table for indexing the patch position
%        edgeMap: input low-res edge
%   Output:     
%        unary: unary potential
%        edgePots: pairwise potential
%        edgeEnds: edge index
% 
% (c)2014 Jun Xie


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

% matrix transform from c++ conversion to matlab conversion
[B,I] = sort(edgeEnds(:, 1));
tmp = edgeEnds(I, :);
edgeEnds(:, 1) = tmp(:, 2);
edgeEnds(:, 2) = tmp(:, 1);

edgePots = edgePots(I, :, :);

edgePots = reshape(edgePots, [size(edgePots,1), labels, labels]);
