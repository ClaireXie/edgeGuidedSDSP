
function runEvaluation(inputFile, scale, highres, scaleFact, border, print2File)

% Evaluation Fuction
%   Input: 
%        inputFile:     input file name
%        scale:         upsampling scales
%        highres:       upsampled result
%        scaleFact:     scaling factor for converting from depth to disparity
%                       (only for middlebury dataset)
%        border:        border to exclude in the evalutation
%        print2File:    indicate whether to print the resultto file
% 
% (c)2016 Jun Xie

para.K=[0.01 0.02];
para.win=fspecial('gaussian', 11, 1.5);

image_file = ['inputs/', inputFile, '_clean.png'];
if exist(image_file,'file')
    gtFile=image_file;
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

gt = double(gt)/scaleFact;
output = highres/scaleFact;

para.l = max(max(gt))-min(min(gt));

rmseV = calc_rmse_roi(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border));

ssimV = ssim(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border),para.K,para.win,para.l);

percentErrorV = calc_err(output(border+1:end-border,border+1:end-border), ... 
    gt(border+1:end-border,border+1:end-border), 1);

if print2File
    fid = fopen(sprintf('outputs/%s_result.txt', inputFile), 'w');
    fprintf(fid, [inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf(fid, [inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf(fid, [inputFile,': percErr= ', num2str(percentErrorV),'\n']);
    fclose(fid); 
else
    fprintf([inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf([inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf([inputFile,': percErr= ', num2str(percentErrorV),'\n']);
end


