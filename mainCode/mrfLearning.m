
function [highres,edges] = mrfLearning(name, indexn, w1, w2, localSize, scale, threshold, show)

% MRF learning with smooth constraint
%   Input: 
%        name: input file name
%        indexn: index of the file in the namelist (for indexing the GT)
%        w1: unary weight -- weights for the shock edge map
%        w2: pairwise weight
%        localSize: overlap size between patches
%        scale: upsampling scales
%        threshold: threshold for exacting the edge (using canny edge)
%        show: indicate whether to visualize the results
%   Output:     
%        highres: upsampled result
%        edges: constructed high-res edges
% 
% (c)2014 Jun Xie

addpath('mainCode/');
addpath('funcs/');

dictName = sprintf('dictionaries/patchData_%d_high', scale);
load(dictName);
addpath(genpath('utils/'));
addpath('mainCode/');
addpath('funcs/');

psize = sqrt(size(lowdata,2));
psize_high = psize;
half = (psize+1)/2;

inputFile = name{indexn};

% for middlebury data
if indexn <= 4
    original = imread(['inputs/', inputFile, '_clean.png']);
    %crop the original for downsampling
    sz = size(original);
    sz = sz - mod(sz, scale);
    original = original(1:sz(1), 1:sz(2));

    input = imresize(original,1/scale,'nearest');
    
% for laser data
else
    load (['inputs/', inputFile]);
    if indexn == 5
        input = bilateralOMA(D);
    else
        input = D;
    end
end

input = double(input);

% nn interpolation
low = imresize(input,scale,'nearest');
edgesl = edge_2010(low,'canny',threshold);

% shock filter parameters
para.dt = 0.1;
para.h = 1;
para.iter = 20;   % iter = 120 for scale = 8
para.lam = 0.00;
para.lam_tld = 1;
para.a = 0.4;
para.theta = pi/1000;
para.smooth = 0;

useMex = 1;
tic;

low0 = real(shock(low,para.iter,para.dt,para.h,'cmp',[para.lam,para.lam_tld,para.a])); 
edgesl0 = edge_2010(low0,'canny',0.1);

[candidateH, candidateHTrans, index, diff] = ...
    genCandidate (edgesl, edgesl0, highdataU, lowdataTrans, highdataTrans, ... 
    psize, localSize, w1);
fprintf('[generate candidates] '); toc; tic;
disp('---------------------------------');

[unary edgePots edgeEnds] = createGraphMex(candidateH, candidateHTrans, index, diff, edgesl);

fprintf('[construct graph] '); toc; tic;
disp('---------------------------------');

[labelling, energy1, energy2] = infer(w2, unary, edgePots, edgeEnds);
fprintf('[infer] '); toc; tic;
disp('---------------------------------');


%reconstruct
output = zeros(size(edgesl));
count = ones(size(edgesl));
for i = 1:length(index)
    ii = (index(i,1))-half+1;
    jj = (index(i,2))-half+1;
    output(ii:ii+psize_high-1,jj:jj+psize_high-1) = ...
        output(ii:ii+psize_high-1,jj:jj+psize_high-1)+...
        reshape(candidateH(i, labelling(i),:),psize_high, psize_high);
    count(ii:ii+psize_high-1,jj:jj+psize_high-1) = ...
        count(ii:ii+psize_high-1,jj:jj+psize_high-1)+1;
end
fprintf('[reconstruct] '); toc; tic;
disp('---------------------------------');

% averaging and thresholding
output = output./count;
output = (output>0.2);


if (show)
    figure;imshow(uint8(output*255));
end

edges = output;

% jopint bilateral super resolution
fprintf('super-resolution...\n');

offset = -1;   % very important
if (useMex)
    % super -resolution (mex code, simplified)
    highres = blup_lowPathMex(edges, input, scale, offset);
else 
    % super -resolution (original code, more accurate)
    highres = blup_lowPath(scale, input, output);
end

fprintf('[bilteral filtering] '); toc; tic;
disp('---------------------------------');

% save the result
if indexn <= 4  
    if (show)
        figure;imshow(uint8(highres));
        imwrite(uint8(highres),['outputs/', inputFile, '2_', num2str(scale), '.png']);
    end
else
    
    if (~exist('outputs','dir'))
        mkdir('outputs');
    end
    save(['outputs/',inputFile, '_SRout.mat'],'highres');
    s = highres(11:end-11,11:end-11);
    tmp = (highres-min(s(:)))/(max(s(:))-min(s(:)));
    if (show)
        figure;imshow(uint8(tmp*255));
        imwrite(uint8(tmp*255),['outputs/',inputFile, '_SR.png']);
    end
    imwrite(uint8(output*255),['outputs/', inputFile, '2_edge_', num2str(scale), '.png']);
end









