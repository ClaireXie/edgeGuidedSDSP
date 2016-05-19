
function [candidateH, candidateHTrans, index, diff1] = ...
           genCandidate (edgeMap, edgeMapShock, highData, lowdataTrans, ... 
           highdataTrans, psize, localSize, w1, useANN)
       
% Function for generating candidate patches for dataset
%   Input: 
%        edgeMap: input low-res edge
%        edfeMapShock: input low-res edge after shock filtering
%        highData: high-res patches from dataset
%        lowdataTrans: low-res patches (distance transformed) from dataset
%        highdataTrans: high-res patches (distance transformed) from dataset
%        psize: window size
%        localSize: overlap size between patches
%        w1: unary weight -- weights for the shock edge map
%        useANN: a boolean var to indicate if the efficient KNN
%        implementation is enabled
%   Output:     
%        candidateH: high-res candidate patches
%        candidateHTrans: high-res candidate patches(distance transformed)
%        index: lookup table for indexing the patch position
%        diff: intensity difference between the candidate and the input
% 
% (c)2014 Jun Xie

half = (psize+1)/2;
numCandidates = 5;

edgeslTmp = edgeMap;
index=[];
num=1;

% extract patches
fprintf('extract patches from the input...\n');
for i = half+1:size(edgeMap,1)-half
    for j = half+1:size(edgeMap,2)-half
        
        if (edgeslTmp(i,j) == 1)
            
            % original patches
            patch = edgeMap(i-half+1:i+half-1, j-half+1:j+half-1);
            % shock patches
            patch0 = edgeMapShock(i-half+1:i+half-1, j-half+1:j+half-1);
                        
            % perform distance transform
            patchTrans = bwdist(patch);
            
            % get the query data for knn
            query(num,:) = double([patchTrans(:); w1*patch0(:)]);
            
                    
            ii = i; jj = j;
            index(num,:) = [ii,jj];

            
            % avoid repeating too many neighbor pixels
            edgeslTmp(i-localSize:i+localSize, j-localSize:j+localSize) = ...
                zeros(2*localSize+1);
            
            num = num+1;
            
            % for visualization only
            %{
            output(i-half+1:i+half-1, j-half+1) = 1;
            output(i-half+1:i+half-1, j+half-1) = 1;
            output(i-half+1, j-half+1:j+half-1) = 1;
            output(i+half-1, j-half+1:j+half-1) = 1;
            %}
            
        end
    end
end


% search for candidates, do knn search
fprintf('finding the knn...\n');
sdata = [lowdataTrans w1*highData];


% add patch match

%matlab implementation of knn
if (~useANN)
    [idx, diff] = knnsearch(sdata, query, 'K', numCandidates);
else
% currently implemented with a more efficient kd tree representation 
% with ANN wrapper
    [idx, diff] = knn_mex(sdata, query, numCandidates);
end

diff1 = [];
% collect candidate patches
for j =0   %:4
    for i = 1:size(idx,1)
        candidateH(i,j*numCandidates+1:(j+1)*numCandidates,:) = ... 
            shift_patch(highData(idx(i,:), :), j, psize);
        candidateHTrans(i,j*numCandidates+1:(j+1)*numCandidates,:) = ... 
            shift_patch(highdataTrans(idx(i,:), :), j, psize);
    end
    diff1 = [diff1 diff];
end


function output = shift_patch(input, dir, w)

if (dir == 0)
    output = input;
    return;
end

input = reshape(input, [], w, w);
output = zeros(size(input));

% move up 
if (dir == 1)
   output(:, 1:w-1, :) =  input(:, 2:w, :);
% move down
elseif (dir == 2)
   output(:, 2:w, :) =  input(:, 1:w-1, :);
% move left
elseif (dir == 3)
   output(:, :, 1:w-1) = input(:, :, 2:w);
elseif (dir == 4) 
   output(:, :, 2:w) = input(:, :, 1:w-1);
end

output = reshape(output, [], w*w);
   
