clear all
close all

nSeats = 6;
nRows = 40;
nNodes = nSeats*nRows;
nStates = 2;
adj = zeros(nNodes);
for r = 1:nRows
	for s = 1:nSeats
		if s < nSeats
			adj(s + nSeats*(r-1),(s+1) + nSeats*(r-1)) = 1; % Passenger behind
		end
		if r < nRows
			adj(s + nSeats*(r-1),s + nSeats*(r)) = 1; % Passenger to the right
		end
	end
end
adj = adj+adj'; % Symmetrize
edgeStruct = UGM_makeEdgeStruct(adj,nStates);

if 0
% Drawing the graph
s = 1;
for r = 1:nRows
	fprintf('(%3d)-(%3d)-(%3d)-(%3d)-(%3d)-(%3d)\n',s,s+1,s+2,s+3,s+4,s+5);
	s = s+nSeats;
	
	if r < nRows
		fprintf('  |     |     |     |     |     |\n');
	end
end
end

alpha = .5;
beta = 2;
nodePot = [ones(nNodes,1) alpha*ones(nNodes,1)];
edgePot = repmat([beta 1;1 beta],[1 1 edgeStruct.nEdges]);

optimalDecoding = UGM_Decode_Junction(nodePot,edgePot,edgeStruct)
[nodeBel,edgeBel,logZ] = UGM_Infer_Junction(nodePot,edgePot,edgeStruct);
nodeBel

samples = UGM_Sample_Junction(nodePot,edgePot,edgeStruct);
figure;
imagesc(reshape(samples(:,1)',nSeats,nRows)');
figure;
imagesc(reshape(samples(:,2)',nSeats,nRows)');
figure;
imagesc(reshape(samples(:,3)',nSeats,nRows)');
figure;
imagesc(reshape(samples(:,4)',nSeats,nRows)');
