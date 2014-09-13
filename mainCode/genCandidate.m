function [candidateH, candidateHTrans, index, diff, output] = ...
           genCandidate (edgeMap, edfeMapShock, highData, lowdataTrans, highdataTrans, scale, psize, localSize, w1)

psize_high = psize;
half = (psize+1)/2;
numCandidates = 5;

edgeslTmp = edgeMap;
index=[];
num=1;

output=edgeMap;
% extract patches
fprintf('extract patches from the input...\n');
for i = half+1:size(edgeMap,1)-half
    for j = half+1:size(edgeMap,2)-half
        
        if (edgeslTmp(i,j) == 1)
            
            % original patches
            patch = edgeMap(i-half+1:i+half-1, j-half+1:j+half-1);
            % shock patches
            patch0 = edfeMapShock(i-half+1:i+half-1, j-half+1:j+half-1);
            
            %{
            [edgelist, labelededgeim] = edgelink(patch, 5);
            
            if (size(edgelist,2)>1)
                k=labelededgeim(half,half);
                [edgelist0, labelededgeim0] = edgelink(patch0, 5);
                patch=(labelededgeim==k);
                patch0=(labelededgeim0==k);
            end
            %}
            
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

% for visualization purpose
%figure;imshow(uint8(output*255));

% search for candidates, do knn search
fprintf('finding the knn...\n');
sdata = [lowdataTrans w1*highData];
[idx, diff] = knnsearch(sdata, query, 'K', numCandidates);

% contruct graph structure
% TODO put this piece of code to mex function
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



