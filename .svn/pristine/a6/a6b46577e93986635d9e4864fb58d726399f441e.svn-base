% constructing the edge map via a linear learnt map result/coeffs;

clc;clear;close all;

addpath('mainCode\');
addpath('funcs\');
addpath('edges\');

load result/patchData_4_high;
load result/coeffs;

scale=4;
psize=sqrt(size(lowdata,2));
%psize_high=psize*scale;
psize_high=psize;
half=(psize+1)/2;
w1=5;
localSize=1;

inputFile='cones';
original=imread(['inputs/', inputFile, '_clean.png']);
%crop the original for downsampling
sz = size(original);
sz = sz - mod(sz, scale);
original = original(1:sz(1), 1:sz(2));

input=imresize(original,1/scale,'nearest');
%input=imread('inputs/teddy_tiny.png');
input=double(input);
low=imresize(input,scale,'nearest');
edgesl = edge(low,'canny',0.1);

low0=imresize(input,scale,'bicubic');


para.dt=0.1;
para.h=1;
para.iter=20;   %120 for scale=8
para.lam=0.00;
para.lam_tld=1;
para.a=0.4;
para.theta=pi/1000;
para.smooth=0;

low0=real(shock(low,para.iter,para.dt,para.h,'cmp',[para.lam,para.lam_tld,para.a])); 
%low0=real(shock(low0,para.iter,para.dt,para.h,'cmp',[para.lam,para.lam_tld,para.a])); 
edgesl0 = edge(low0,'canny',0.1);


edgeslTmp=edgesl;
num=1;
output=zeros(size(input)*scale);
count=ones(size(input)*scale);


%linear learning

for i=half+1:size(edgesl,1)-half
    for j=half+1:size(edgesl,2)-half
        if (edgeslTmp(i,j)==1)
            patch=edgesl(i-half+1:i+half-1, j-half+1:j+half-1);
            patch_predict=patch(:)'*w;
            ii=i-half+1;
            jj=j-half+1;
            output(ii:ii+psize_high-1,jj:jj+psize_high-1)=...
                output(ii:ii+psize_high-1,jj:jj+psize_high-1)+...
            reshape(patch_predict,psize_high, psize_high);
    
            count(ii:ii+psize_high-1,jj:jj+psize_high-1)=...
                count(ii:ii+psize_high-1,jj:jj+psize_high-1)+1;  
            
            % avoid repeating too many neighbor pixels
            edgeslTmp(i-localSize:i+localSize, j-localSize:j+localSize)=0;
        end
    end
end

output=output./count;
output0=(output>0.2);
figure;imshow(uint8(output*255));
figure;imshow(uint8(output0*255));

imwrite(uint8(output*500),'teddy_linearLear.png');


% super-resolution
%{
fprintf('super-resolution...\n');
highres=blup_low(scale, input, output);
figure;imshow(uint8(highres));
imwrite(uint8(highres),['outputs/', inputFile, '_', num2str(scale), '.png']);
%}



