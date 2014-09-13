function template=computeDistTemp(n)

A=zeros(n^2);
% 8 connected neighbor
for i=1:n^2
    if (mod(i,n)~=0)
        A(i,i+1)=1;
    end
    if (i<n^2-n)
         A(i,i+n)=1;
    end
    if (mod(i,n)~=1 && i<n^2-n)
        A(i,i+n-1)=sqrt(2);
    end
    if (mod(i,n)~=0 && i<n^2-n)
        A(i,i+n+1)=sqrt(2);
    end
end


A=A+transpose(A);

C=A;
A=min(A,1);
center=(n^2+1)/2;
[costs,paths] = dijkstra(A,C,center);
%[costs,paths] = dijkstra(A,C);

template=zeros(n,n,n^2);
NHOOD=[0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);


for i=1:n^2
    p=paths{1,i};
    % do not need to count the starting point and ending point
    for j=2:length(p)-1
        x=ceil(p(j)/n);
        y=mod(p(j)-1, n)+1;
        template(x,y,i)=1;
    end
    % make 4 neighbor connected

    template(:,:,i) = imdilate(template(:,:,i), se);
    
    x1=ceil(p(1)/n);
    y1=mod(p(1)-1, n)+1;
    template(x1,y1,i) = 1;
    
    x2=ceil(p(end)/n);
    y2=mod(p(end)-1, n)+1;
    template(x2,y2,i) = 1;
end





