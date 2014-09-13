function [A, C] = constructGraphEdge(n)

% TODO: make A and C sparse
% A: edge in the graph
% C: edge weifht in the graph

A = zeros(n^2);
C = A;
for i = 1:n^2    
    if (mod(i,n) ~= 0)
        A(i,i+1) = 1;
        C(i,i+1) = 1;
    end
    if (i < n^2-n)
         A(i,i+n) = 1;
         C(i,i+n) = 1;
    end
    if (mod(i,n) ~= 1 && i < n^2-n)
        A(i,i+n-1) = 1;
        C(i,i+n-1) = sqrt(2);
    end
    if (mod(i,n) ~= 0 && i < n^2-n)
        A(i,i+n+1) = 1;
        C(i,i+n+1) = sqrt(2);
    end
end