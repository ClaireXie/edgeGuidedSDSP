clear all
close all

%% Test 1: Independent Nodes
nodePot = [1 3
    9 1
    1 3
    9 1]

nSamples = 100;
for s = 1:nSamples
   samples(1,s) = sampleDiscrete(nodePot(1,:)/sum(nodePot(1,:))); 
   samples(2,s) = sampleDiscrete(nodePot(2,:)/sum(nodePot(2,:))); 
   samples(3,s) = sampleDiscrete(nodePot(3,:)/sum(nodePot(3,:))); 
   samples(4,s) = sampleDiscrete(nodePot(4,:)/sum(nodePot(4,:))); 
end
figure;
imagesc(samples');
ylabel('Question');
xlabel('Student');
title('Test 1');

%% Test 2: Dependent nodes

nNodes = 4;
nStates = 2;

% Make adjacency matrix
adj = zeros(nNodes,nNodes);
adj(1,2) = 1;
adj(2,1) = 1;
adj(2,3) = 1;
adj(3,2) = 1;
adj(3,4) = 1;
adj(4,3) = 1;
fprintf('This is what the adjacency matrix looks like:\n');
adj

% Make edgeStruct
fprintf('This is what the edgeStruct looks like:\n');
edgeStruct = UGM_makeEdgeStruct(adj,nStates)

fprintf('Here is the number of states vector:\n');
edgeStruct.nStates

fprintf('Here are the edges:\n');
edgeEnds = edgeStruct.edgeEnds

fprintf('Here are the nodes associated with edge 2:\n');
edge = 2;
edgeStruct.edgeEnds(edge,:)

fprintf('Here are the edge numbers connected to node 3:\n');
n = 3;
edges = UGM_getEdges(n,edgeStruct);

fprintf('Here are the neighbours of node %d\n',n);
nodes = edgeStruct.edgeEnds(edges,:);
nodes(nodes ~= n)

% Make node potentials
fprintf('Here are the node potentials\n');
nodePot = [1 3
    9 1
    1 3
    9 1]

% Make edge 
maxState = max(edgeStruct.nStates);
edgePot = zeros(maxState,maxState,edgeStruct.nEdges);
for e = 1:edgeStruct.nEdges
   edgePot(:,:,e) = [2 1 ; 1 2]; 
end
fprintf('Here are the edge potentials\n');
edgePot

fprintf('Computing Optimal Decoding...\n');
optimalDecoding = UGM_Decode_Exact(nodePot,edgePot,edgeStruct)

fprintf('Computing Node Marginals, Edge Marginals, and Log of Normalizing Constant...\n');
[nodeBel,edgeBel,logZ] = UGM_Infer_Exact(nodePot,edgePot,edgeStruct);
nodeBel
fprintf('Normalization constant: %f\n',exp(logZ));

edgeStruct.maxIter = 100;
fprintf('Generating %d samples from the model...\n',edgeStruct.maxIter);
samples = UGM_Sample_Exact(nodePot,edgePot,edgeStruct);

figure;
imagesc(samples')
ylabel('Question');
xlabel('Student');
title('Test 2');
