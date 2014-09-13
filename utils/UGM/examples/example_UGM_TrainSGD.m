%% Load Noisy X
clear all
close all

load X.mat
[nRows,nCols] = size(X);
nNodes = nRows*nCols;
nStates = 2;
nInstances = 100;

% Make 100 noisy X instances
y = int32(1+X);
y = reshape(y,[1 1 nNodes]);
y = repmat(y,[nInstances 1 1]);

X = reshape(X,1,1,nNodes);
X = repmat(X,[nInstances 1 1]);
X = X + randn(size(y))/2;

figure;
for i = 1:4
subplot(2,2,i);
imagesc(reshape(X(i,1,:),nRows,nCols));
colormap gray
end
suptitle('Examples of Noisy Xs');

%% Make edgeStruct

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
nEdges = edgeStruct.nEdges;

%% Make Xnode, Xedge, infoStruct, initialize weights

% Add bias and Standardize Columns
tied = 1;
Xnode = [ones(nInstances,1,nNodes) UGM_standardizeCols(X,tied)];
nNodeFeatures = size(Xnode,2);

% Make nodeMap
nodeMap = zeros(nNodes,nStates,nNodeFeatures,'int32');
for f = 1:nNodeFeatures
    nodeMap(:,1,f) = f;
end

% Make Xedge
sharedFeatures = [1 0];
Xedge = UGM_makeEdgeFeatures(Xnode,edgeStruct.edgeEnds,sharedFeatures);
nEdgeFeatures = size(Xedge,2);

% Make edgeMap
f = max(nodeMap(:));
edgeMap = zeros(nStates,nStates,nEdges,nEdgeFeatures,'int32');
for edgeFeat = 1:nEdgeFeatures
   edgeMap(1,1,:,edgeFeat) = f+edgeFeat;
   edgeMap(2,2,:,edgeFeat) = f+edgeFeat;
end

%% Evaluate with random parameters

figure;
nParams = max([nodeMap(:);edgeMap(:)]);
w = randn(nParams,1);
for i = 1:4
    subplot(2,2,i);
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    nodeBel = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
    imagesc(reshape(nodeBel(:,2),nRows,nCols));
    colormap gray
end
suptitle('Loopy BP node marginals with random parameters');
fprintf('(paused)\n');
pause

%% Train with Loopy Belief Propagation for 3 iterations

maxIter = 3; % Number of passes through the data set

w = zeros(nParams,1);
options.maxFunEvals = maxIter;
funObj = @(w)UGM_CRF_NLL(w,Xnode,Xedge,y,nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP);
w = minFunc(funObj,w,options);

figure;
for i = 1:4
    subplot(2,2,i);
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    nodeBel = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
    imagesc(reshape(nodeBel(:,2),nRows,nCols));
    colormap gray
end
suptitle('Loopy BP node marginals with truncated minFunc parameters');
fprintf('(paused)\n');
pause

%% Train with Stochastic gradient descent for the same amount of time
stepSize = 1e-4;
w = zeros(nParams,1);
for iter = 1:maxIter*nInstances
    i = ceil(rand*nInstances);
    funObj = @(w)UGM_CRF_NLL(w,Xnode(i,:,:),Xedge(i,:,:),y(i,:),nodeMap,edgeMap,edgeStruct,@UGM_Infer_LBP);
    [f,g] = funObj(w);
    
    fprintf('Iter = %d of %d (fsub = %f)\n',iter,maxIter*nInstances,f);
    
    w = w - stepSize*g;
end

figure;
for i = 1:4
    subplot(2,2,i);
    [nodePot,edgePot] = UGM_CRF_makePotentials(w,Xnode,Xedge,nodeMap,edgeMap,edgeStruct,i);
    nodeBel = UGM_Infer_LBP(nodePot,edgePot,edgeStruct);
    imagesc(reshape(nodeBel(:,2),nRows,nCols));
    colormap gray
end
suptitle('Loopy BP node marginals with truncated  stochastic gradient parmaeters');
fprintf('(paused)\n');
pause


