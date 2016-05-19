function output_name = extractPatches(D, scaling, psize, num_rot)

% Script for extracting the training set of patches based on self
% similarity
%   Input: 
%        D: Input low-res image
%        scaling: upsampling scales
%        psize: patch size
%        num_rot: number of rotation samples
%   Output:     
%        output name of the dictionary
% 
% (c)2016 Jun Xie


disp('========================')
fprintf('Collect high-res and low-res edge pairs for Self-similarity...\n');
 
psize_half = (psize+1)/2;

preAlloc = zeros(10000, psize*psize);
allocationTime = 0;
batchSize = 5000;
highdata = preAlloc;
lowdata = preAlloc;
highdataTrans = preAlloc;
lowdataTrans = preAlloc;
showEdges = 0;

%-----------Training (Extracting Edges from the training set) ---------%
% crop the data
[r,c] = size(D);
D = D(1:end-rem(r,scaling),1:end-rem(c,scaling));
low = imresize(D,1/scaling,'nearest');

% downsample
low = imresize(low,scaling,'nearest');

edgesh = edge(double(D),'canny',0.03);
edgesl = edge(double(low),'canny',0.03);

% extract patches
% generate a mask for only edge region
mask = edgesl;
cIt = 1;
for i = 1:size(edgesl,1)-psize
    for j = 1:size(edgesl,2)-psize
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
                % random sample one out of three dirs
                if (num_rot > 1)
                    for r = 2:num_rot
                        theta = 180/num_rot*(r-1);
                        rpatch_high = extractRotatedPatch(edgesh, [i+psize_half-1, j+psize_half-1], psize, psize, theta);
                        rpatch_low = extractRotatedPatch(edgesl, [i+psize_half-1, j+psize_half-1], psize, psize, theta);
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
                end


                if (showEdges)
                    % plot the result
                    subplot(1,2,1); imshow(reshape(lowdata(cIt-1,:),21,21));
                    subplot(1,2,2); imshow(reshape(highdata(cIt-1,:),21,21));
                    keyboard;
                end

            end

        end  
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

% save to (temp) result
output_name = sprintf('dictionaries/patchData_%d_high_temp', scaling);

save(output_name, ... 
    'highdataU','lowdataU', 'highdataTrans', 'lowdataTrans');


