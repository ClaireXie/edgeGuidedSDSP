function [nodeBel, edgeBel, logZ] = UGM_Infer_Junction(nodePot,edgePot,edgeStruct,ordering)

debug = 0;

[nNodes,maxState] = size(nodePot);
nEdges = size(edgePot,3);
edgeEnds = edgeStruct.edgeEnds;
V = edgeStruct.V;
E = edgeStruct.E;
nStates = edgeStruct.nStates;

if nargin < 4
    ordering = 1:nNodes;
end

%% Triangulate graph and find cliques
adj = UGM_VE2adj(V,E,edgeEnds);
adj_new = adj;
cliques = cell(0,1);
for n = ordering(:)'

    % Add extra
    neighbors = find(adj(:,n));
    for n1 = neighbors(:)'
        for n2 = neighbors(:)'
            if n1 ~= n2
                adj(n1,n2) = 1;
                adj(n2,n1) = 1;
                adj_new(n1,n2) = 1;
                adj_new(n1,n2) = 1;
            end
        end
    end

    % Add to set of cliques
    clique(n,n) = 1;
    clique(n,neighbors) = 1;

    % Remove node
    adj(n,:) = 0;
    adj(:,n) = 0;

    % If not a subset of an existing clique, add to set of cliques
    cand = [n;neighbors];
    add = 1;
    for c = 1:length(cliques)
        if isempty(setdiff(cand,cliques{c}))
            add = 0;
            break;
        end
    end
    if add
        cliques{end+1,1} = cand;
    end

    if debug
        clf;
        drawGraph(adj_new)
        pause
    end
end
nCliques = length(cliques);

if debug;
    fprintf('Cliques:\n')
    cliques{:}
end

%% Make initial potentials
nodePotMissing = ones(nNodes,1);
edgePotMissing = ones(nEdges,1);
for c = 1:nCliques
    fprintf('Clique %d: %s\n',c,sprintf(' %d',cliques{c}));
    nodes = cliques{c};

    % Initialize clique potentials
    cliquePot{c} = ones([nStates(nodes)' 1]);

    ind = cell(length(nodes),1);
    for nodeInd = 1:length(nodes)
        ind{nodeInd} = 1:nStates(nodes(nodeInd));
    end

    for nodeInd = 1:length(nodes)
        n = nodes(nodeInd);
        if nodePotMissing(n)
            fprintf('Including node potential for node %d in clique %d\n',n,c);
            nodePotMissing(n) = 0;

            ind_sub = ind;
            for s = 1:nStates(n)
                ind_sub{nodeInd} = s;
                cliquePot{c}(ind_sub{:}) = cliquePot{c}(ind_sub{:})*nodePot(n,s);
            end
        end

        edges = E(V(n):V(n+1)-1);
        for e = edges(:)'
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);
            if ismember(n1,nodes) && ismember(n2,nodes) && edgePotMissing(e)
                fprintf('Including edge potential for edge %d-%d in clique %d\n',edgeEnds(e,1),edgeEnds(e,2),c);
                edgePotMissing(e) = 0;

                ind_sub = ind;
                nodeInd1 = find(nodes==n1);
                nodeInd2 = find(nodes==n2);
                for s1 = 1:nStates(n1)
                    for s2 = 1:nStates(n2)
                        ind_sub{nodeInd1} = s1;
                        ind_sub{nodeInd2} = s2;
                        cliquePot{c}(ind_sub{:}) = cliquePot{c}(ind_sub{:})*edgePot(s1,s2,e);
                    end
                end
            end
        end
    end
    fprintf('\n');
end
if debug
    fprintf('Clique Potentials:\n');
    cliquePot{:}
end

%% Makes edges between cliques that satisfy RIP
weights = zeros(0,3);
for c1 = 1:nCliques
    for c2 = c1+1:nCliques
        sep = length(intersect(cliques{c1},cliques{c2}));
        if sep > 0
            weights(end+1,:) = [c1 c2 -sep];
        end
    end
end
edges = find(minSpan(nCliques,weights));
cliqueEdges = weights(edges,1:2);
nCliqueEdges = size(cliqueEdges,1);
[V,E] = UGM_makeEdgeVE(cliqueEdges,nCliques);
for e = 1:nCliqueEdges
    separators{e,1} = intersect(cliques{cliqueEdges(e,1)},cliques{cliqueEdges(e,2)});
end
if debug
    cliqueEdges
    fprintf('Separators:\n');
    separators{:}
end

%% Message Passing

% Count number of neighbors
nNeighbors = zeros(nCliques,1);
for e = 1:nCliqueEdges
    nNeighbors(cliqueEdges(e,1)) = nNeighbors(cliqueEdges(e,1))+1;
    nNeighbors(cliqueEdges(e,2)) = nNeighbors(cliqueEdges(e,2))+1;
end

% Add all leafs to initial queue
Q = find(nNeighbors == 1);

sent = zeros(nCliqueEdges*2,1);
waiting = ones(nCliqueEdges*2,1);
messages = cell(nCliqueEdges*2,1);
while ~isempty(Q)
    c = Q(1);
    Q = Q(2:end);

    wait = waiting(V(c):V(c+1)-1);
    sending = sent(V(c):V(c+1)-1);

    nWaiting = sum(wait==1);

    if nWaiting == 0
        % Send final messages
        for sendEdge = [double(V(c))+find(sending==0)-1]'
            sent(sendEdge) = 1;
            [messages,waiting,nei] = send(c,sendEdge,cliques,cliquePot,messages,waiting,nStates,cliqueEdges,V,E,separators);
            if nNeighbors(nei) == 1 || nNeighbors(nei) == 0
                Q = [Q;nei];
            end
        end
    else
        remainingEdge = V(c)+find(wait==1)-1;
        sent(remainingEdge) = 1;
        [messages,waiting,nei] = send(c,remainingEdge,cliques,cliquePot,messages,waiting,nStates,cliqueEdges,V,E,separators);
        nNeighbors(nei) = nNeighbors(nei)-1;
        if nNeighbors(nei) == 1 || nNeighbors(nei) == 0
            Q = [Q;nei];
        end
    end
end
%messages{:}

%% Compute cliqueBel
cliqueBel = cell(nCliques,1);
for c = 1:nCliques
    nodes = cliques{c};
    ind = cell(length(nodes),1);
    for nodeInd = 1:length(nodes)
        ind{nodeInd} = 1:nStates(nodes(nodeInd));
    end

    % Multiply cliquePot by all incoming messages
    cb = cliquePot{c};
    edges = E(V(c):V(c+1)-1);
    for e = edges(:)'
        if c == cliqueEdges(e,2)
            msg = messages{e};
        else
            msg = messages{e+nCliqueEdges};
        end


        ind_sub = ind;
        sepLength = length(separators{e});
        sep = zeros(sepLength,1);
        s = cell(length(sep),1);
        for n = 1:sepLength
            s{n,1} = 1;
            sep(n) = find(nodes==separators{e}(n));
        end
        while 1
            for nodeInd = 1:length(sep)
                ind_sub{sep(nodeInd)} = s{nodeInd};
            end
            cb(ind_sub{:}) = cb(ind_sub{:})*msg(s{:});

            for nodeInd = 1:sepLength
                s{nodeInd} = s{nodeInd} + 1;
                if s{nodeInd} <= nStates(sep(nodeInd))
                    break;
                else
                    s{nodeInd} = 1;
                end
            end
            if nodeInd == length(sep) && s{end} == 1
                break;
            end
        end
    end
    cb = cb./sum(cb(:));
    cliqueBel{c} = cb;
end
%cliqueBel{:}

%% Compute nodeBel
nodeBel = zeros(size(nodePot));
nodeBelMissing = ones(nNodes,1);
for c = 1:nCliques
    cb = cliqueBel{c};

    nodes = cliques{c};
    ind = cell(length(nodes),1);
    for nodeInd = 1:length(nodes)
        ind{nodeInd} = 1:nStates(nodes(nodeInd));
    end

    for nodeInd = 1:length(nodes)
        n = nodes(nodeInd);
        if nodeBelMissing(n)
            nodeBelMissing(n) = 0;

            ind_sub = ind;
            for s = 1:nStates(n)
                ind_sub{nodeInd} = s;
                slice = cb(ind_sub{:});
                nodeBel(n,s) = sum(slice(:));
            end
        end
    end
end

%% Compute edgeBel
if nargout > 1
    edgeBel = zeros(size(edgePot));
    edgeBelMissing = ones(nEdges,1);
    for c = 1:nCliques
        cb = cliqueBel{c};

        nodes = cliques{c};
        ind = cell(length(nodes),1);
        for nodeInd = 1:length(nodes)
            ind{nodeInd} = 1:nStates(nodes(nodeInd));
        end

        for e = 1:nEdges
            n1 = edgeEnds(e,1);
            n2 = edgeEnds(e,2);
            n1Ind = find(n1==nodes);
            n2Ind = find(n2==nodes);
            if edgeBelMissing(e) && ~isempty(n1Ind) && ~isempty(n2Ind)
                edgeBelMissing(e) = 0;

                ind_sub = ind;
                for s1 = 1:nStates(n1)
                    for s2 = 1:nStates(n2)
                        ind_sub{n1Ind} = s1;
                        ind_sub{n2Ind} = s2;
                        slice = cb(ind_sub{:});
                        edgeBel(s1,s2,e) = sum(slice(:));
                    end
                end
            end
        end
    end
end

if nargout > 2
   %% Compute sepBel
   sepBel = cell(nCliqueEdges,1);
   for e = 1:nCliqueEdges
       c = cliqueEdges(e,1);
       cb = cliqueBel{c};
       nodes = cliques{c};
       ind = cell(length(nodes),1);
       for nodeInd = 1:length(nodes)
           ind{nodeInd} = 1:nStates(nodes(nodeInd));
       end
       
        sepLength = length(separators{e});
        sep = zeros(sepLength,1);
        s = cell(length(sep),1);
        for n = 1:sepLength
            s{n,1} = 1;
            sep(n) = find(nodes==separators{e}(n));
        end
        sb = ones([nStates(separators{e})' 1]);
        while 1
            for nodeInd = 1:length(sep)
                ind{sep(nodeInd)} = s{nodeInd};
            end
            sb(s{:}) = sum(cb(ind{:}));

            for nodeInd = 1:sepLength
                s{nodeInd} = s{nodeInd} + 1;
                if s{nodeInd} <= nStates(sep(nodeInd))
                    break;
                else
                    s{nodeInd} = 1;
                end
            end
            if nodeInd == length(sep) && s{end} == 1
                break;
            end
        end
        sepBel{e} = sb./sum(sb(:));
   end
   
   % Node Energy
   Energy1 = 0;
   for n = 1:nNodes
       Energy1 = Energy1 - sum(nodeBel(n,1:nStates(n)).*log(nodePot(n,1:nStates(n))));
   end
   % Edge Energy
   Energy2 = 0;
   for e = 1:nEdges
      n1 = edgeEnds(e,1);
      n2 = edgeEnds(e,2);
      eb = edgeBel(1:nStates(n1),1:nStates(n2),e);
      ep = edgePot(1:nStates(n1),1:nStates(n2),e);
      Energy2 = Energy2 - sum(eb(:).*log(ep(:)));
   end
   % Separator Entropy
   Entropy1 = 0;
   for e = 1:nCliqueEdges
      Entropy1 = Entropy1 + sum(sepBel{e}(:).*log(sepBel{e}(:)));
   end
   % Clique Entropy
   Entropy2 = 0;
   for c = 1:nCliques
       Entropy2 = Entropy2 - sum(cliqueBel{c}(:).*log(cliqueBel{c}(:)));
   end
   F = (Energy1+Energy2) - (Entropy1+Entropy2);
   logZ = -F;
end


end

%% Message passing function
function [messages,waiting,nei] = send(c,e,cliques,cliquePot,messages,waiting,nStates,edgeEnds,V,E,separators)
nEdges = size(edgeEnds,1);
edge = E(e);
if c == edgeEnds(edge,1)
    nei = edgeEnds(edge,2);
else
    nei = edgeEnds(edge,1);
end
fprintf('Sending from %d to %d\n',c,nei);

% Opposite edge is no longer waiting
for tmp = V(nei):V(nei+1)-1
    if tmp ~= e && E(tmp) == E(e)
        waiting(tmp) = 0;
    end
end

e = edge;

nodes = cliques{c};
for nodeInd = 1:length(nodes)
    ind{nodeInd} = 1:nStates(nodes(nodeInd));
end

% Compute Product of clique potential with all incoming messages except
% along e
temp = cliquePot{c};
neighbors = E(V(c):V(c+1)-1);
for e2 = neighbors(:)'
    if e ~= e2
        ind_sub = ind;
        sepLength = length(separators{e2});
        sep = zeros(sepLength,1);
        s = cell(length(sep),1);
        for n = 1:sepLength
            s{n,1} = 1;
            sep(n) = find(nodes==separators{e2}(n));
        end
        while 1
            for nodeInd = 1:length(sep)
                ind_sub{sep(nodeInd)} = s{nodeInd};
            end
            if c == edgeEnds(e2,2)
                temp(ind_sub{:}) = temp(ind_sub{:})*messages{e2}(s{:});
            else
                temp(ind_sub{:}) = temp(ind_sub{:})*messages{e2+nEdges}(s{:});
            end

            for nodeInd = 1:length(sep)
                s{nodeInd} = s{nodeInd} + 1;
                if s{nodeInd} <= nStates(sep(nodeInd))
                    break;
                else
                    s{nodeInd} = 1;
                end
            end
            if nodeInd == length(sep) && s{end} == 1
                break;
            end
        end
    end
end

% Sum out over all variables except separator set
sepLength = length(separators{e});
sep = zeros(sepLength,1);
s = cell(length(sep),1);
for n = 1:sepLength
    s{n,1} = 1;
    sep(n) = find(nodes==separators{e}(n));
end
newm = ones([nStates(sep)' 1]);
ind_sub = ind;
while 1
    for nodeInd = 1:length(sep)
        ind_sub{sep(nodeInd)} = s{nodeInd};
    end
    slice = temp(ind_sub{:});
    newm(s{:}) = sum(slice(:));

    for nodeInd = 1:length(sep)
        s{nodeInd} = s{nodeInd} + 1;
        if s{nodeInd} <= nStates(sep(nodeInd))
            break;
        else
            s{nodeInd} = 1;
        end
    end
    if nodeInd == length(sep) && s{end} == 1
        break;
    end
end
newm = newm./sum(newm(:));

if c == edgeEnds(e,2)
    messages{e+nEdges} = newm;
else
    messages{e} = newm;
end
end