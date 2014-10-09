
function [candidateH, candidateHTrans, index, diff] = ...
           genCandidate (edgeMap, edfeMapShock, highData, lowdataTrans, ... 
           highdataTrans, psize, localSize, w1)
       
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
            patch0 = edfeMapShock(i-half+1:i+half-1, j-half+1:j+half-1);
                        
            % perform distance transform
            patchTrans = bwdist(patch);
            
            % get the query data for knn
            query(num,:) = double([patchTrans(:); w1*patch0(:)]);
            
            ii = i; jj = j;

            index(num,:) = [ii,jj];
            patchLow(num,:) = patch(:)';
            patchHigh(num,:) = patch0(:)';
            
            % avoid repeating too many neighbor pixels
            edgeslTmp(i-localSize:i+localSize, j-localSize:j+localSize) = ...
                zeros(2*localSize+1);
            
            num = num+1;
            
            % for visualization
            output(i-half+1:i+half-1, j-half+1) = 1;
            output(i-half+1:i+half-1, j+half-1) = 1;
            output(i-half+1, j-half+1:j+half-1) = 1;
            output(i+half-1, j-half+1:j+half-1) = 1;
            
        end
    end
end


% search for candidates, do knn search
fprintf('finding the knn...\n');
sdata = [lowdataTrans w1*highData];
[idx, diff] = knnsearch(sdata, query, 'K', numCandidates);

% contruct graph structure
% older implementation, currently implemented in mex
%{
fprintf('contruct graph structure and collect candidate patches...\n');
structure = zeros(size(index,1));
structure = sparse(structure);
for i = 1:size(index,1)
    % determine whether two patches have overlaps
    x0 = index(i,1);
    y0 = index(i,2);
    for j = i+1:size(index,1)
        xx = index(j,1);
        yy = index(j,2);
        
        if (abs(x0-xx) < 21 && abs(y0-yy) < 21)
            structure(i,j) = 1;
        end
    end
end
toc; tic;
%}

% collect candidate patches
for i = 1:size(idx,1)
    candidateH(i,:,:) = highData(idx(i,:), :);
    candidateHTrans(i,:,:) = highdataTrans(idx(i,:), :);
end



