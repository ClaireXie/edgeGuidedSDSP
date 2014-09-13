function [ unary edgePots edgeEnds] = createGraph(candidateH, candidateHTrans, index, diff, structure, edgeMap, half)

% define the unary term
nNodes = size(index,1);
unary = diff.^2;

c = size(candidateH,2);
psize = sqrt(size(candidateH,3));

[rows, cols, s] = find(structure);

fprintf('construct the edge weights...\n');

for i = 1:size(cols,1)
    a = rows(i);
    b = cols(i);
    
    tmpMap = zeros(size(edgeMap));
    tmpMap(index(a,1)-half+1:index(a,1)+half-1, index(a,2)-half+1:index(a,2)+half-1) = 1;
    tmpMap(index(b,1)-half+1:index(b,1)+half-1, index(b,2)-half+1:index(b,2)+half-1) = ...
        tmpMap(index(b,1)-half+1:index(b,1)+half-1, index(b,2)-half+1:index(b,2)+half-1)+1;
    mask = (tmpMap == 2);
    maskA = mask(index(a,1)-half+1:index(a,1)+half-1, index(a,2)-half+1:index(a,2)+half-1);
    maskB = mask(index(b,1)-half+1:index(b,1)+half-1, index(b,2)-half+1:index(b,2)+half-1);
    
    [ra,ca] = find(maskA);
    [rb,cb] = find(maskB);
    
    current = reshape(candidateHTrans(a,:,:),c, psize,psize);
    %current=reshape(candidateH(a,:,:),c, psize,psize);
    
    current = current(:,min(ra):max(ra),min(ca):max(ca));
    
    for k=1:c 
        
        next = reshape(candidateHTrans(b,k,:),psize,psize);
        %next=reshape(candidateH(b,k,:),psize,psize);
        
        next = next(min(rb):max(rb),min(cb):max(cb));
        
        overlapDiff = reshape(current,c,[])-repmat(next(:)', c, 1);
        edgePots(i, :, k) = sum(overlapDiff.^2, 2);
    end
    
    edgeEnds(i, :) = [a b];
        
end
