% mrf learning with smooth constraint

function [highres,edges] = mrfLearning2(name, indexn, w1, w2, localSize, scale, threshold, show)

addpath('mainCode/');
addpath('funcs/');
addpath('edges/');

load dictionaries/patchData_4_high;
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
edgesl = edge(low,'canny',threshold);  

% bicubic interpolation
% meaningless (just for comparison)
low0 = imresize(input,scale,'bicubic');

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
edgesl0 = edge(low0,'canny',0.1);


[candidateH, candidateHTrans, index, diff, mask] = ...
    genCandidate (edgesl, edgesl0, highdataU, lowdataTrans, highdataTrans, scale, psize, localSize, w1);
fprintf('[generate candidates] '); toc; tic;
disp('---------------------------------');

[unary edgePots edgeEnds] = createGraphMex(candidateH, candidateHTrans, index, diff, edgesl);

%save('result/graph_5.mat','unary', 'edgePots', 'edgeEnds');
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
% for temporally usage
edges = output;


% save the temp result for debugging
save('tempResult.mat', 'edges', 'input', 'scale');

fprintf('super-resolution...\n');

%{
%crop the original for downsampling
colorImg = imread(['inputs/', inputFile, '_color.png']);
sz = [size(colorImg, 1), size(colorImg, 2)];
sz = sz - mod(sz, scale);
colorImg = colorImg(1:sz(1), 1:sz(2), :);
%colorImg = imresize(colorImg, sz);
edges = edge(rgb2gray(colorImg),'canny',0.2); 

addpath(genpath('sketchTokens/'));
opts=struct('nPos',100,'nNeg',80,'modelFnm','modelSmall','nTrees',20);
model=stTrain(opts);
addpath(genpath('toolbox/'));
st = stDetect( colorImg, model ); 
E = stToEdges( st, 1 );
edges = E>0.75;
rmpath(genpath('toolbox/'));
cd .. 
rmpath(genpath('sketchTokens/'));

tic;
%}

offset = -1;    %-1
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
    save(['outputs/',inputFile, '_SRout.mat'],'highres');
    s = highres(11:end-11,11:end-11);
    tmp = (highres-min(s(:)))/(max(s(:))-min(s(:)));
    if (show)
        figure;imshow(uint8(tmp*255));
        imwrite(uint8(tmp*255),['outputs/',inputFile, '_SR.png']);
    end
    imwrite(uint8(output*255),['outputs/', inputFile, '2_edge_', num2str(scale), '.png']);
end









