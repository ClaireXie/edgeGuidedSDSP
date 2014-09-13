function template=findPath(A, C, ending, n)

center = (n^2+1)/2;
[costs,p] = dijkstra(A,C,center,ending);

template = zeros(n,n);

NHOOD = [0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);

%draw the path
for j = 2:length(p)-1
    x = ceil(p(j)/n);
    y = mod(p(j)-1, n)+1;
    template(x,y) = 1;
end
% make 4 neighbor connected

template = imdilate(template, se);
template((n+1)/2,(n+1)/2) = 1;

x2 = ceil(p(end)/n);
y2 = mod(p(end)-1, n)+1;
template(x2,y2) = 1;
