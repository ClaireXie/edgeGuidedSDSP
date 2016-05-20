% Script for extracting the training set of patches
% Training script (Extracting Edges from the training set)
% (c)2016 Jun Xie

clc; clear; close all;
addpath('funcs/');

disp('========================')
fprintf('collect high-res and low-res edge pairs...\n');

scaling = 4;
psize = 21; 
psize_half = (psize+1)/2;
fnames = dir('800/*.mat');

preAlloc = zeros(200000, psize*psize);
allocationTime = 0;
batchSize = 50000;
highdata = preAlloc;
lowdata = preAlloc;
highdataTrans = preAlloc;
lowdataTrans = preAlloc;
iter_prev = 0;

showEdges = 0;

%-----------Training (Extracting Edges from the training set) ---------%
cIt = 1;
for k = 1:length(fnames)
    D = [];
    name = fnames(k).name;
    load(['800/', name]);
    
    % extract edges for high/low-res
    % crop the data
    [r,c] = size(D);
    D = D(1:end-rem(r,scaling),1:end-rem(c,scaling));
    low = imresize(D,1/scaling,'nearest');
    low = imresize(low,scaling,'nearest');
    
    edgesh = edge(double(D),'canny',0.03);
    edgesl = edge(double(low),'canny',0.03);
    
    % extract patches
    % generate a mask for only edge region
    mask = edgesl;
    temp = edgesh;
    delta = ((temp-edgesl)>0);

    for i = 1:2:size(edgesl,1)-psize
        for j = 1:2:size(edgesl,2)-psize
            if (mask(i+psize_half-1,j+psize_half-1) ~= 0)
                plow = edgesl(i:i+psize-1,j:j+psize-1);
                phigh = edgesh(i:i+psize-1,j:j+psize-1);
                
                plowtrans = bwdist(plow);
                phightrans = bwdist(phigh);

                % determine whether the edge maps are consistent
                t1 = sum(plow(:));
                t2 = sum(phigh(:));
                if (t1>2*scaling && t2>2*scaling && abs(t2-t1)<0.8*psize)
                    
                    highdata(cIt,:)= phigh(:)';
                    lowdata(cIt,:) = plow(:)';
                    highdataTrans(cIt,:)= phightrans(:)';
                    lowdataTrans(cIt,:) = plowtrans(:)';

                    cIt = cIt+1;
                    if (cIt > size(preAlloc, 1)+batchSize*allocationTime - 3) 
                        % if pre-allocation is not enogh, double the size
                        highdata = [highdata; zeros(batchSize, psize*psize)];
                        lowdata = [lowdata; zeros(batchSize, psize*psize)];
                        highdataTrans = [highdataTrans; zeros(batchSize, psize*psize)];
                        lowdataTrans = [lowdataTrans; zeros(batchSize, psize*psize)];
                        allocationTime = allocationTime + 1;
                    end
                    
                    %also added reotated patches
                    for r = 2:3
                        theta = 60*(r-1);
                        rpatch_high = extractRotatedPatch(edgesh, [j+psize_half-1, i+psize_half-1], psize, psize, theta);
                        rpatch_low = extractRotatedPatch(edgesl, [j+psize_half-1, i+psize_half-1], psize, psize, theta);
                        rpatch_high = rpatch_high>0;
                        rpatch_low = rpatch_low>0;
                        
                        highdata(cIt,:)= rpatch_high(:)';
                        lowdata(cIt,:) = rpatch_low(:)';

                        phightrans = bwdist(rpatch_high);
                        plowtrans = bwdist(rpatch_low);
                        highdataTrans(cIt,:)= phightrans(:)';
                        lowdataTrans(cIt,:) = plowtrans(:)';
                        
                        cIt=cIt+1;
                    end
                    
                       
                    if (showEdges)
                        % plot the result
                        subplot(1,2,1); imshow(reshape(lowdata(cIt,:),21,21));
                        subplot(1,2,2); imshow(reshape(highdata(cIt,:),21,21));
                    end
                end
            end  
        end
    end
    
    iter=floor(k/length(fnames)*10);
    if iter~=iter_prev
        fprintf('%d%% ',iter*10);
        iter_prev=iter;
        
        [~, ia, ~] = unique([highdata(1:cIt-1, :), lowdata(1:cIt-1, :)],'rows');
        highdata(1:length(ia), :) = highdata(ia,:);
        lowdata(1:length(ia), :) = lowdata(ia,:);
        highdataTrans(1:length(ia), :) = highdataTrans(ia,:);
        lowdataTrans(1:length(ia), :) = lowdataTrans(ia,:);

        highdata(length(ia)+1:end, :) = 0;
        lowdata(length(ia)+1:end, :) = 0;
        highdataTrans(length(ia)+1:end, :) = 0;
        lowdataTrans(length(ia)+1:end, :) = 0; 

        cIt = length(ia)+1;
    end
end

%de-allocate
highdata(cIt:end, :) = [];
lowdata(cIt:end, :) = [];
highdataTrans(cIt:end, :) = [];
lowdataTrans(cIt:end, :) = [];

fprintf('\n');

[C, ia, ic] = unique([highdata, lowdata],'rows');
highdataU = highdata(ia,:);
lowdataU = lowdata(ia,:);
highdataTrans = highdataTrans(ia,:);
lowdataTrans = lowdataTrans(ia,:);
fprintf('saving the patch data...\n');

% save the result
save(sprintf('dictionaries/patchData_%d_high', scaling), ... 
    'highdataU','lowdataU', 'highdataTrans', 'lowdataTrans');
