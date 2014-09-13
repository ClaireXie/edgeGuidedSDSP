clc;clear;close all;
scaling=2;
psize=11;
psize_half=(psize+1)/2;
fnames = dir('800/*.mat');

highdata=[];
lowdata=[];
psize_high=psize*scaling;
iter_prev=0;

for k=1:length(fnames)
    D=[];
    name=fnames(k).name;
    load(['800/', name]);
    
    %extract edges for high/low-res
    % crop the data
    [r,c]=size(D);
    D=D(1:end-rem(r,scaling),1:end-rem(c,scaling));
    low=imresize(D,1/scaling,'nearest');
    
    edgesh = edge(double(D),'canny',0.03);
    edgesl = edge(double(low),'canny',0.08);
    
    % extract patches
    
    % generate a mask for only edge region
    mask=edgesl;
    %mask=imdilate(edgesl, strel('disk',(psize+1)/2));
    temp=imresize(edgesh, 1/scaling, 'nearest');
    delta=((temp-edgesl)>0);

    cIt = size(highdata,1)+1;
    for i=1:size(edgesl,1)-psize
        for j=1:size(edgesl,2)-psize
            if (mask(i+psize_half-1,j+psize_half-1)~=0)
                ii = (i-1) * scaling + 1;
                jj = (j-1) * scaling + 1;
                plow=edgesl(i:i+psize-1,j:j+psize-1);
                phigh=edgesh(ii:ii+psize_high-1,jj:jj+psize_high-1);

                % determine whether the edge maps are consistent
                if (sum(plow(:))>2 && sum(phigh(:))>2*scaling)
                    if sum(sum(delta(i:i+psize-1,j:j+psize-1)))==0
                        highdata(cIt,:)= phigh(:)';
                        lowdata(cIt,:) = plow(:)';
                        cIt=cIt+1;
                    end  
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
highdataU=highdata(ia,:);
lowdataU=lowdata(ia,:);

fprintf('saving the patch data...\n');
save('patchData_4','highdata','lowdata','highdataU','lowdataU');
