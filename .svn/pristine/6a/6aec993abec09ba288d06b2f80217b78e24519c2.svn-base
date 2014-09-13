%clc; close all; clear;
function output=blup_low(scale, input, edge0)


depthLow=input;

%---------------------------------------%
depthMedium=imresize(depthLow, scale, 'bicubic');
%edgesh = edge(double(depthMedium),'canny');
%edgesh = edge(double(depthOrg),'canny');
%edgesh=imresize(edge0, scale, 'bicubic');
edgesh=edge0;
%---------------------------------------%

%either declaire as global or assign value individually
global window sigma_d sigma_c;
%window=5;
%sigma_d=0.5;
window_half=floor(window/2);
offset=1; % very important
depthHigh=zeros(size(depthMedium));

%generate rastered line template
template=computeDistTemp(scale*2*window_half+1);

for i=scale*window_half+1 : size(edgesh,1)-scale*window_half-offset
    for j=scale*window_half+1 : size(edgesh,2)-scale*window_half-offset
        s=0;normalize=0;
               
         avg=0;
         candidate=[];
         gdist=[];
                
        for ii=i-scale*window_half:i+scale*window_half
            for jj=j-scale*window_half: j+scale*window_half
                supportE=edgesh(i-scale*window_half:i+scale*window_half,...
                    j-scale*window_half:j+scale*window_half);
                
                x=ii-(i-scale*window_half)+1;
                y=jj-(j-scale*window_half)+1;
               
                
                if (mod(ii+offset,scale)==0 && mod(jj+offset,scale)==0)
                    index=y+(x-1)*(scale*2*window_half+1); 
                    ind = find(template(:,:,index)~=0);
                    
                    if sum(supportE(ind))==0 || (edgesh(i,j)==1 && sum(supportE(ind))==1)
                        
                        avg=avg+depthLow((ii+offset)/scale,(jj+offset)/scale);
                        candidate=[candidate, depthLow((ii+offset)/scale,(jj+offset)/scale)];
                        d=sqrt((i/scale-ii/scale)^2+(j/scale-jj/scale)^2); 
                        g_d=normpdf(d, 0, sigma_d);
                        gdist=[gdist, g_d];
                        
                        %s=s+depthLow((ii+offset)/scale,(jj+offset)/scale)*g;
                        %normalize=normalize+g;
                    end       
                end
            end
        end
        
        %if (normalize~=0)
        if (~isempty(candidate))
             avg=avg/length(candidate);
            for canNum=1:length(candidate)
                g_c=normpdf(candidate(canNum), avg, sigma_c);
                s=s+candidate(canNum)*g_c*gdist(canNum);
                normalize=normalize+g_c*gdist(canNum);
            end
            depthHigh(i,j)=s/normalize;
        else
            depthHigh(i,j)=depthMedium(i,j);
        end
    end
end

output=depthHigh;
%figure;imshow(uint8(depthHigh));
%figure;imshow(uint8(depthHigh+edgesh*255));

