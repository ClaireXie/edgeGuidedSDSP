function output = blup_lowPath(scale, input, edge0)


depthLow = input;

%---------------------------------------%
depthMedium = imresize(depthLow, scale, 'bicubic');
%edgesh = edge(double(depthMedium),'canny');
%edgesh = edge(double(depthOrg),'canny');
%edgesh=imresize(edge0, scale, 'bicubic');
edgesh = edge0;
%---------------------------------------%

% TODO: need to avoid the global vars
global window sigma_d;
window_half = floor(window/2);
%----------------------------------------%
offset = 1; % very important
%----------------------------------------%
depthHigh = zeros(size(depthMedium));

%generate rastered line template
% find the shortest path from one point to all other pixels 
% so template0 is of size n*n*(n*n)
template0 = computeDistTemp(scale*2*window_half+1);
[structA, structC0] = constructGraphEdge(scale*2*window_half+1);
[row, col] = find(structA == 1);
structA = structA+transpose(structA);
structC0 = structC0+transpose(structC0);

% processing the edges (dilation)
% To avoid disconnection between edge points
NHOOD = [0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);
edgeMap = edge0+imdilate(edge0, se)+1;


% TODO:
% 1. fix the problem as stated in the loop
% 2. adpapt this function into mex function
for idx = 1:size(edgesh,1)*size(edgesh,2)
        
        % be careful about the order
        j = ceil(idx/size(edgesh,1));
        i = mod(idx-1, size(edgesh,1))+1;
       
        % don't count the pixel on the border
        if (i<scale*window_half+1 || i>size(edgesh,1)-scale*window_half-offset ||...
                j<scale*window_half+1 || j>size(edgesh,2)-scale*window_half-offset)
            continue;
        end
         
        %{
        % be careful about the order
        idx=(j-1)*size(edgesh,1)+i;
        %}
        
        % accumulated vars
        s = 0; 
        normalize = 0;
        
        edgeTmp = edgeMap(i-scale*window_half:i+scale*window_half, ...
            j-scale*window_half:j+scale*window_half);        
        centerLoc = edgeMap(idx);
        
        % if centerLoc == 2, then this pixel is not originally an 
        % edge point (it is from dilation--> along the edge)
        
        if (centerLoc == 2)
            % compute the edge weights
            structC = constructGraphWeight(row, col, edgeTmp, scale*2*window_half+1, structC0);
        end
               
        % bilateal filtering
        for ii = i-scale*window_half:i+scale*window_half
            for jj = j-scale*window_half: j+scale*window_half
                
                supportE = edgesh(i-scale*window_half:i+scale*window_half,...
                    j-scale*window_half:j+scale*window_half);
                
                x = ii-(i-scale*window_half)+1;
                y = jj-(j-scale*window_half)+1;
                
                % only compute pixel values when it has values in the 
                % low-res image
                if (mod(ii+offset,scale) == 0 && mod(jj+offset,scale) == 0)
                    index = y+(x-1)*(scale*2*window_half+1); 
                       
                    % if this patch contain edges
                    % only close to the edge pixels
                    
                    % TODO: 
                    % i feel it might have some problems
                    % if q is along the edge, we still need to recompute
                    % the edge weights?
                    
                    % another problem:
                    % it could be possible that both p and q are not on/along the
                    % edge, but the shortst path are along the edge
                    %----------------------------------------%
                    if (centerLoc == 2)
                        template = findPath(structA, structC, index, scale*2*window_half+1);
                    else
                        % otherwise, just use the pre-computed template
                        % matrix directly
                        template = template0(:,:,index);
                    end
                    %----------------------------------------%
                   
                    ind = find(template ~= 0);
                    
                    if sum(supportE(ind)) == 0 || (edgesh(i,j) == 1 && sum(supportE(ind)) == 1)
                         d = sqrt((i/scale-ii/scale)^2+(j/scale-jj/scale)^2);
                         g = normpdf(d, 0, sigma_d);
                         s = s+depthLow((ii+offset)/scale,(jj+offset)/scale)*g;
                         normalize = normalize+g;
                    end       
                end
            end
        end
        
        
        if (normalize ~= 0)
            depthHigh(idx) = s/normalize;
        else
            depthHigh(idx) = depthMedium(i,j);
        end
    %end
end

output = depthHigh;


