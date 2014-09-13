function runEvaluation(inputFile, scale, highres, edges, scaleFact, border, threshold, print2File)

global fid;
offset=0;
para.K=[0.01 0.02];
para.win=fspecial('gaussian', 11, 1.5);

gtFile=['inputs/', inputFile, '_clean.png'];
outputFile=['outputs/v2/s', num2str(scale), '/', inputFile, '2_', num2str(scale), '.png'];


gt=imread(gtFile);
%crop the original for downsampling
sz = size(gt);
sz = sz - mod(sz, scale);
gt = gt(1:sz(1), 1:sz(2));

edgesGt = edge(gt,'canny',threshold);   %0.1
gt=double(gt)/scaleFact;
output=highres/scaleFact;
output(1+offset:end,1+offset:end)=output(1:end-offset, 1:end-offset);

para.l=max(max(gt))-min(min(gt));



rmseV=calc_rmse(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border));

ssimV=ssim(output(border+1:end-border,border+1:end-border),...
    gt(border+1:end-border,border+1:end-border),para.K,para.win,para.l);

edgeError=sum(sum(abs(edgesGt(border+1:end-border,border+1:end-border)-...
    edges(border+1:end-border,border+1:end-border))))...
    /((size(edgesGt,1)-2*border)*(size(edgesGt,2)-2*border));

if print2File
    fprintf(fid, [inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf(fid, [inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf(fid, [inputFile,': edgeErr= ', num2str(edgeError),'\n']);
else
    fprintf([inputFile,': RMSE=    ', num2str(rmseV),'\n']);
    fprintf([inputFile,': SSIM=    ', num2str(ssimV),'\n']);
    fprintf([inputFile,': edgeErr= ', num2str(edgeError),'\n']);
end
