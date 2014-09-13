function [labelling, energy1, energy2]=infer(weight, unary, edgePots0, edgeEnds)

% optimization setting
maxIter = 100;
useMex = 1;

% graph
nNodes = size(unary,1);
nState = size(unary,2);
edgePots = weight*edgePots0;
[V,E] = UGM_makeEdgeVE(edgeEnds,nNodes,useMex);
edgeStruct.edgeEnds = int32(edgeEnds);
edgeStruct.V = V;
edgeStruct.E = E;
edgeStruct.nNodes = nNodes;
edgeStruct.nEdges = size(edgeEnds,1);
edgeStruct.nStates = int32(repmat(nState,nNodes,1));
edgeStruct.useMex = useMex;
edgeStruct.maxIter = maxIter;

edgePots = permute(edgePots,[2, 3, 1]);
edgePots = exp(-1*edgePots/100);
nodePots = exp(-1*unary/100);

%do we need normalization??
%{
for i=1:size(edgePots,3)
    edgePots(:,:,i)=edgePots(:,:,i)./repmat(sum(edgePots(:,:,i),2),1,10);
end
nodePots=nodePots./repmat(sum(nodePots,2),1,10);
%}

fprintf('[Inference] Decoding using L-BP...\n');

% decoding
edgeStruct.nStates = int32(edgeStruct.nStates);
LBPDecoding = UGM_Decode_LBP(nodePots,edgePots,edgeStruct);
labelling = LBPDecoding;

% output energy
energy2 = 0;
for i = 1:size(unary,1)
    energy2 = energy2+unary(i,labelling(i));
end
for i = 1:size(edgePots0,1)
    ii = edgeEnds(i, 1);
    jj = edgeEnds(i, 2);
    energy2 = energy2+weight*edgePots0(i,labelling(ii), labelling(jj));
end


energy1 = 0;
for i = 1:size(unary,1)
    energy1 = energy1+unary(i,1);
end
for i = 1:size(edgePots0,1)
    energy1 = energy1+weight*edgePots0(i, 1, 1);
end

