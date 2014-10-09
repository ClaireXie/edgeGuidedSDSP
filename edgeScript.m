
% script for extracting the training set of patches
% (c)2014 Jun Xie

clc; clear; close all;

disp('========================')
fprintf('collection high-res and low-res edge pairs...\n');

% set the parameters below
%--------------------------------%
scaling = 4;
psize = 21; 
psize_half = (psize+1)/2;
fnames = dir('800/*.mat');

highdata = [];
lowdata = [];
highdataTrans = [];
lowdataTrans = [];
iter_prev = 0;

showEdges = 0;

%-----------Training (Extracting Edges from the training set) ---------%

for k = 1:length(fnames)
    D = [];
    name = fnames(k).name;
    load(['800/', name]);
    
    % extract edges for high/low-res
    % crop the data
    [r,c] = size(D);
    D = D(1:end-rem(r,scaling),1:end-rem(c,scaling));
    low = imresize(D,1/scaling,'nearest');
    
    % new added
    low = imresize(low,scaling,'nearest');
    
    edgesh = edge(double(D),'canny',0.03);
    edgesl = edge(double(low),'canny',0.03);
    
    % extract patches
    % generate a mask for only edge region
    mask = edgesl;
    temp = edgesh;
    delta = ((temp-edgesl)>0);

    cIt = size(highdata,1)+1;
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
                        lowdataTrans(cIt,:) = phightrans(:)'; 
                    
                    if (showEdges)
                        % plot the result
                        subplot(1,2,1); imshow(reshape(lowdata(cIt,:),21,21));
                        subplot(1,2,2); imshow(reshape(highdata(cIt,:),21,21));
                    end
                    
                    cIt=cIt+1;
                end

            end  
        end
    end
    
    iter=floor(k/length(fnames)*10);
    if iter~=iter_prev
        fprintf('%d%% ',iter*10);
        iter_prev=iter;
    end
    
end

fprintf('\n');

temp=[highdata, lowdata];
[C, ia, ic] = unique(temp,'rows');
highdataU = highdata(ia,:);
lowdataU = lowdata(ia,:);
highdataTrans = highdataTrans(ia,:);
lowdataTrans = lowdataTrans(ia,:);
fprintf('saving the patch data...\n');

% save the result
save('result/patchData_4_high','highdata','lowdata','highdataU','lowdataU', 'highdataTrans', 'lowdataTrans');
