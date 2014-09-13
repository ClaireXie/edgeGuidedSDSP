clc; close all; clear;

window=5;
sigma_d=0.5;
window_half=floor(window/2);
scale=4;

%generate gaussian table
gtable=zeros(window,window);
for i=1:window
    for j=1:window
        d=sqrt((i-window_half)^2+(j-window_half)^2);
        gtable(i,j)=normpdf(d, 0, sigma_d);
    end
end
%normalize=sum(sum(gtable));

depthLow=imread('inputs/disp2_tiny.png');
depthLow=double(depthLow);
depthMedium=imresize(depthLow, [375,450], 'bicubic');
depthHigh=depthMedium;

depthOrg=imread('inputs/cones_clean.png');
%depthMedium2=imresize(depthLow, [375,450], 'bicubic');
edgesh = edge(double(depthOrg),'canny');

%generate rastered line template (only for 5 by 5)
tmp(:,:,1)=[1 1 0; 
            1 1 1; 
            0 1 1];
        
tmp(:,:,2)=[0 1 0; 
            0 1 1; 
            0 0 1];

tmp(:,:,3)=[0 0 1; 
            0 0 1; 
            0 0 1];

tmp(:,:,4)=[0 0 0; 
            1 1 0; 
            0 1 1];
        
tmp(:,:,5)=[0 0 0; 
            0 1 1; 
            0 1 1];
        
tmp(:,:,6)=[0 0 0; 
            0 0 1; 
            0 0 1];
        
tmp(:,:,7)=[0 0 0; 
            0 0 0; 
            1 1 1];
        
tmp(:,:,8)=[0 0 0; 
            0 0 0; 
            0 1 1];

template=zeros(window, window, window*window);
template(1:3,1:3,1:3)=tmp(:,:,1:3);
template(1:3,3:5,4)=fliplr(tmp(:,:,2));
template(1:3,3:5,5)=fliplr(tmp(:,:,1));

template(1:3,1:3,6:8)=tmp(:,:,4:6);
template(1:3,3:5,9)=fliplr(tmp(:,:,5));
template(1:3,3:5,10)=fliplr(tmp(:,:,4));

template(1:3,1:3,11:12)=tmp(:,:,7:8);
template(1:3,1:3,13)=zeros(3, 3);
template(3,3,13)=1;
template(1:3,3:5,14)=fliplr(tmp(:,:,8));
template(1:3,3:5,15)=fliplr(tmp(:,:,7));

template(3:5,1:3,16)=flipud(tmp(:,:,4));
template(3:5,1:3,17)=flipud(tmp(:,:,5));
template(3:5,1:3,18)=flipud(tmp(:,:,6));
template(3:5,3:5,19)=fliplr(flipud(tmp(:,:,5)));
template(3:5,3:5,20)=fliplr(flipud(tmp(:,:,4)));

template(3:5,1:3,21)=flipud(tmp(:,:,1));
template(3:5,1:3,22)=flipud(tmp(:,:,2));
template(3:5,1:3,23)=flipud(tmp(:,:,3));
template(3:5,3:5,24)=fliplr(flipud(tmp(:,:,2)));
template(3:5,3:5,25)=fliplr(flipud(tmp(:,:,1)));


for i=window_half+1:size(depthMedium,1)-window_half
    for j=window_half+1:size(depthMedium,2)-window_half
        support=depthMedium(i-window_half:i+window_half,j-window_half:j+window_half);
        supportE=edgesh(i-window_half:i+window_half,j-window_half:j+window_half);
        s=0;normalize=0;
        for ii=1:window
            for jj=1:window
                index=jj+(ii-1)*window;           
                ind = find(template(:,:,index)~=0);
                if sum(supportE(ind))==0 
                    s=s+support(ii,jj)*gtable(ii,jj);
                    normalize=normalize+gtable(ii,jj);
                % if the point is on the edge, use blending mode
                elseif (supportE(ii, jj)==1 && ii==window_half+1 && jj==window_half+1)
                    s=s+sum(sum(support(ii-1:ii+1, jj-1:jj+1)))*gtable(ii,jj)/9;
                    normalize=normalize+gtable(ii,jj);
                end
            end
        end
        %s=sum(sum(support.*gtable));
        depthHigh(i,j)=s/normalize;
    end
end

figure;imshow(uint8(depthMedium));
figure;imshow(uint8(depthHigh));
figure;imshow(uint8(depthHigh+edgesh*255));

