
function highRes = blup_lowPathMex(edges, input, scale, offset)
% Mexed function for bilateral upsampling
%   Input: 
%        edges: input low-res edge
%        input: low-res image
%        scale: upsampling scales
%        offset: edge offset (fixed)
%   Output:     
%        highRes: high-res reconstructed depth image
% 
% (c)2014 Jun Xie

global window sigma_d;

depthMedium = imresize(input, scale, 'nearest');

window_half = floor(window/2);
template = computeDistTemp(scale*2*window_half+1);

[structA, structC0] = constructGraphEdge(scale*2*window_half+1);
structA = structA+transpose(structA);
structC = structC0+transpose(structC0);
G =sparse(structC);
 
NHOOD = [0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);
edgeMap = edges+imdilate(edges, se)+1;

% change from matlab index to c++ index (except template)
edges = edges';
depthMedium = depthMedium';
edgeMap = edgeMap';
input = input';

highRes = bilateralUpsample(input, depthMedium, double(edges), template, ... 
    window, sigma_d, scale, offset, G, nnz(G), edgeMap);

% convert back to matlab index
highRes = (reshape(highRes, size(edges, 1), size(edges, 2)))';
