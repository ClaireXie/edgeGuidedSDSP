function  [y] = UGM_Decode_AlphaExpansionBetaShrink(nodePot, edgePot, edgeStruct, decodeFunc, betaSelect, y)
% INPUT
% nodePot(node,class)
% edgePot(class,class,edge) where e is referenced by V,E (must be the same
% between feature engine and inference engine)
%
% OUTPUT
% nodeLabel(node)

[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;
maxState = max(nStates);

% Initialize
if nargin < 4
	decodeFunc = @UGM_Decode_GraphCut;
end
if nargin < 5
	betaSelect = 0;
end
if nargin < 6
	[junk y] = max(nodePot,[],2);
end
y = int32(y);
if edgeStruct.useMex
	pot = UGM_LogConfigurationPotentialC(int32(y),nodePot,edgePot,int32(edgeEnds));
else
	pot = UGM_LogConfigurationPotential(y,nodePot,edgePot,edgeStruct.edgeEnds);
end

% Do Alpha-Expansions until convergence
while 1
	y_old = y;
	
	for s1 = 1:maxState
		switch betaSelect
			case 0 % Basic alpha-expansion
				y = applyMove(y,s1,s1,nodePot,edgePot,edgeStruct,decodeFunc);
			case 1 % Randomized selection of beta
				ind = [1:s1-1 s1+1:maxState];
				y = applyMove(y,s1,ind(ceil(rand*(maxState-1))),nodePot,edgePot,edgeStruct,decodeFunc);
			case 2 % Iterate over choices of beta
				%for s2 = [1:s1-1 s1+1:maxState]
				for s2 = 1:maxState
					y = applyMove(y,s1,s2,nodePot,edgePot,edgeStruct,decodeFunc);
				end
			case 3 % Iterate over choices of alpha
				%for s2 = [1:s1-1 s1+1:maxState]
				for s2 = 1:maxState
					y = applyMove(y,s2,s1,nodePot,edgePot,edgeStruct,decodeFunc);
				end
			case 4 % Set beta = alpha-1
				y = applyMove(y,s1,max(s1-1,1),nodePot,edgePot,edgeStruct,decodeFunc);
			case 5 % Set beta = alpha+1
				y = applyMove(y,s1,min(s1+1,maxState),nodePot,edgePot,edgeStruct,decodeFunc);
		end
	end
	
	if all(y==y_old)
		break;
	end
	pot_old = pot;
	if edgeStruct.useMex
		pot = UGM_LogConfigurationPotentialC(int32(y),nodePot,edgePot,int32(edgeEnds));
	else
		pot = UGM_LogConfigurationPotential(y,nodePot,edgePot,edgeStruct.edgeEnds);
	end
	fprintf('logPot = %f, changes = %d\n',pot,sum(y~=y_old));
	if pot-pot_old == 0
		break;
	end
end

end

function y = applyMove(y,s1,s2,nodePot,edgePot,edgeStruct,decodeFunc)
[nNodes,maxStates] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;

fprintf('Expanding %d, shrinking with %d\n',s1,s2);

modifiedNP = zeros(nNodes,2);
modifiedEP = zeros(2,2,nEdges);
if edgeStruct.useMex
	UGM_Decode_AlphaExpansionBetaShrinkC(y,int32(s1-1),int32(s2-1),nodePot,edgePot,edgeEnds,modifiedNP,modifiedEP);
else
	for n = 1:nNodes
		if y(n) == s1
			modifiedNP(n,:) = [nodePot(n,s1) nodePot(n,s2)]; % Keep alpha or shrink with beta
		else
			modifiedNP(n,:) = [nodePot(n,s1) nodePot(n,y(n))]; % Keep current state or switch to alpha
		end
	end
	for e = 1:nEdges
		n1 = edgeEnds(e,1);
		n2 = edgeEnds(e,2);
		y1 = y(n1);
		y2 = y(n2);
		if y1 == s1
			y1 = s2;
		end
		if y2 == s1
			y2 = s2;
		end
		modifiedEP(:,:,e) = [edgePot(s1,s1,e) edgePot(s1,y2,e)
			edgePot(y1,s1,e) edgePot(y1,y2,e)];
		
		% Change energy to make sub-modular (does nothing if edges are metric)
		energy = -log(modifiedEP(:,:,e));
		if y(n1) == s1 && y(n2) == s1
			% (alpha,alpha) edge, decrease energy of staying at (alpha,alpha)
			energy(1,1) = min(energy(1,1),energy(1,2)+energy(2,1)-energy(2,2));
		elseif y(n1) == s1
			% (alpha,~alpha) edge, increase energy of changing to (~alpha,alpha)
			energy(2,1) = max(energy(2,1),energy(1,1)+energy(2,2)-energy(1,2));
		elseif y(n2) == s1
			% (~alpha,alpha) edge, increase energy of changing to (alpha,~alpha)
			energy(1,2) = max(energy(1,2),energy(1,1)+energy(2,2)-energy(2,1));
		else
			% (~alpha,~alpha) edge, decrease of energy of staying at (~alpha,~alpha)
			energy(2,2) = min(energy(2,2),energy(1,2)+energy(2,1)-energy(1,1));
		end
		modifiedEP(:,:,e) = exp(-energy);
	end
end
modifiedES = edgeStruct;
modifiedES.nStates(:) = 2;
ytmp = decodeFunc(modifiedNP,modifiedEP,modifiedES);
y(ytmp == 2 & y == s1) = s2;
y(ytmp == 1) = s1;
end