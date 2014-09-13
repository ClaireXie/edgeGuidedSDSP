
clear all
close all
load X.mat
[nRows,nCols] = size(X);

y = 2-X;
ysub = y(1:16,1:16);
ysub(ysub==2) = 3;
y(1:16,1:16) = ysub;
ysub = y(1:16,17:end);
ysub(ysub==2) = 4;
y(1:16,17:end) = ysub;

Xrgb = ones(nRows,nCols,3);
Xrgb(:,:,2) = Xrgb(:,:,2) - (y==2);
Xrgb(:,:,3) = Xrgb(:,:,3) - (y==2);
Xrgb(:,:,1) = Xrgb(:,:,1) - (y==3);
Xrgb(:,:,3) = Xrgb(:,:,3) - (y==3);
Xrgb(:,:,1) = Xrgb(:,:,1) - (y==4);
Xrgb(:,:,2) = Xrgb(:,:,2) - (y==4);

figure;
imagesc(Xrgb);
title('Original X');

Xrgb = Xrgb + randn(size(Xrgb))/2;
if(exist('imshow')~=0)
figure;
imshow(Xrgb);
end
title('Noisy X');

%% Represent denoising as UGM

[nRows,nCols] = size(X);
nNodes = nRows*nCols;
nStates = 4;

adj = sparse(nNodes,nNodes);

% Add Down Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],repmat(nRows,[1 nCols]),1:nCols); % No Down edge for last row
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+1)) = 1;

% Add Right Edges
ind = 1:nNodes;
exclude = sub2ind([nRows nCols],1:nRows,repmat(nCols,[1 nRows])); % No right edge for last column
ind = setdiff(ind,exclude);
adj(sub2ind([nNodes nNodes],ind,ind+nRows)) = 1;

% Add Up/Left Edges
adj = adj+adj';
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

%% Set up learning problem

% Make labels
y = int32(reshape(y,[1 1 nNodes]));

% Make node features
X = zeros(1,3,nNodes);
X(1,1,:) = reshape(Xrgb(:,:,1),1,1,nNodes);
X(1,2,:) = reshape(Xrgb(:,:,2),1,1,nNodes);
X(1,3,:) = reshape(Xrgb(:,:,3),1,1,nNodes);
tied = 1;
Xnode = UGM_standardizeCols(X,tied);
Xnode = [ones(1,1,nNodes) X];
nNodeFeatures = size(Xnode,2);
nodeMap = zeros(nNodes,nStates,nNodeFeatures,'int32');
p = 1;
for f = 1:nNodeFeatures
	for s = 1:nStates
		nodeMap(:,s,f) = p;
		p = p+1;
	end
end

% Make edge features
sharedFeatures = [1 0 0 0];
Xedge = UGM_makeEdgeFeaturesInvAbsDif(Xnode,edgeStruct.edgeEnds,sharedFeatures);
nEdgeFeatures = size(Xedge,2);
edgeMap = zeros(nStates,nStates,edgeStruct.nEdges,nEdgeFeatures,'int32');
for f = 1:nEdgeFeatures
	for s = 1:nStates
		edgeMap(s,s,:,f) = p;
	end
	p = p+1;
end

nParams = max([nodeMap(:);edgeMap(:)]);
w = zeros(nParams,1);

funObj = @(w)UGM_CRF_PseudoNLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct); % Make objective with new Xedge/edgeMap
UB = inf(nParams,1); % No upper bound on parameters
LB = -inf(nParams,1); % No lower bound on node parameters, 
LB(max(nodeMap(:))+1:end) = 0; % edge parameters must be non-negative  

w = minConf_TMP(funObj,w,LB,UB);

[nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct);
pause
% 
% [w,v] = UGM_initWeights(infoStruct);
% wv = [w(:);v(:)];
% UB = inf(size(wv));
% LB = -inf(size(wv));
% LB(end-length(v(:))+1:end) = 0;
% funObj = @(wv)UGM_CRFpseudoLoss(wv,X,Xedge,y,edgeStruct,infoStruct);
% wv = minConf_TMP(funObj,wv,LB,UB);
% [w,v] = UGM_splitWeights(wv,infoStruct);
% nodePot = UGM_makeCRFNodePotentials(X,w,edgeStruct,infoStruct);
% edgePot = UGM_makeCRFEdgePotentials(Xedge,v,edgeStruct,infoStruct);
