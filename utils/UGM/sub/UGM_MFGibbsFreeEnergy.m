
function [F] = MFGibbsFreeEnergy(nodePot,edgePot,nodeBel,nStates,edgeEnds,V,E)

[nNodes,maxState] = size(nodePot);
nEdges = size(edgeEnds,1);

threshold = 1e-10;

U1 = 0;
U2 = 0;
S1 = 0;

for n = 1:nNodes
    % Local Mean-Field Average Energy Term
    b = nodeBel(n,1:nStates(n));
    U1 = U1 + sum(b .* log(nodePot(n,1:nStates(n))));

    % Mean-Field Entropy Term
    b = nodeBel(n,1:nStates(n));
    b(b < threshold) = 1;
    S1 = S1 + sum(b.* log(b));
end

for e = 1:nEdges
    n1 = edgeEnds(e,1);
    n2 = edgeEnds(e,2);

    b_i = repmat(nodeBel(n1,1:nStates(n1))',1,nStates(n2));
    b_j = repmat(nodeBel(n2,1:nStates(n2)),nStates(n1),1);
    pot_ij = edgePot(1:nStates(n1),1:nStates(n2),e);

    % Pairwise Mean-Field Average Energy Term
    U2 = U2 + sum(b_i(:).*b_j(:).*log(pot_ij(:)));
end

F = - U2 - U1 + S1;
end