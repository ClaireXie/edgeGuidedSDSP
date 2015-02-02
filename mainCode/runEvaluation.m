
function runEvaluation(inputFile, scale, highres, edges, scaleFact, border, ... 
    threshold, print2File, indexn)
% Evaluation Fuction
%   Input: 
%        inputFile: input file name
%        scale: upsampling scales
%        highres: upsampled result
%        edges: upsampled edges
%        scaleFact: scaling factor for converting from depth to disparity
%        (only for middlebury dataset)
%        border: border to exclude in the evalutation
%        threshold: threshold for exacting the GT edge (using canny edge)
%        localSize: overlap size between patches
%        threshold: threshold for exacting the edge (using canny edge)
%        print2File: indicate whether to print the resultto file
%        indexn: index of the file in the namelist (for indexing the GT)
% 
% (c)2014 Jun Xie

fid = fopen(sprintf('outputs/%s_result.txt', inputFile), 'w');
offset=0;
para.K=[0.01 0.02];
para.win=fspecial('gaussian', 11, 1.5);

if (indexn <= 4)
    gtFile=['inputs/', inputFile, '_clean.png'];
    gt=imread(gtFile);
else
    gtFile = ['inputs/gt_laser/', inputFile(1:end-5)];
    load (gtFile);
    gt = D;
end
%crop the original for downsampling
sz = size(gt);
sz = sz - mod(sz, scale);
gt = gt(1:sz(1), 1:sz(2));

edgesGt = edge(gt,'canny',threshold);   %0.1
gt = double(gt)/scaleFact;
output = highres/scaleFact;
output(1+offset:end,1+offset:end) = output(1:end-offset, 1:end-offset);

para.l = max(max(gt))-min(min(gt));

rmseV = calc_rmse(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border));

ssimV = ssim(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border),para.K,para.win,para.l);

percentErrorV = calc_err(output(border+1:end-border,border+1:end-border), ... 
    gt(border+1:end-border,border+1:end-border), 1);

if print2File
    fprintf(fid, [inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf(fid, [inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf(fid, [inputFile,': percErr= ', num2str(percentErrorV),'\n']);
else
    fprintf([inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf([inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf([inputFile,': percErr= ', num2str(percentErrorV),'\n']);
end

% compute error map
% comment the following if you want to save the error image

%errorMapDepth = abs(output(border+1:end-border,border+1:end-border) - ... 
%    gt(border+1:end-border,border+1:end-border));

%errorMapEdge = (edges(border+1:end-border,border+1:end-border) - ... 
%    edgesGt(border+1:end-border,border+1:end-border));

%h1 = figure;imagesc(errorMapDepth);colorbar;title('ErrorMap-depth');truesize;
%h2 = figure;imagesc(errorMapEdge);colorbar;title('ErrorMap-edge');truesize;

%print(h1,'-dpng',sprintf('outputs/errorMapDepth_%s.png', inputFile));
%print(h2,'-dpng',sprintf('outputs/errorMapEdge_%s.png', inputFile));

fclose(fid);  
