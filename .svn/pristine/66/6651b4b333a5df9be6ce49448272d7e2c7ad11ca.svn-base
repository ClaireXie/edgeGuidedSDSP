function template=computeDistPath(edges, n)

A=zeros(n^2);
C=zeros(n^2);
% 8 connected neighbor
for i=1:n^2
    indi=floor(i/n)+1;
    indj=rem(i,n);
    if indj==0
        indj=n;
    end
    
    if (mod(i,n)~=0)
        A(i,i+1)=1;
        C(i,i+1)=max(edges(indi,indj),edges(indi,indj+1));  
    end
    if (i<n^2-n)
         A(i,i+n)=1;
         C(i,i+n)=max(edges(indi,indj),edges(indi+1,indj));
    end
    if (mod(i,n)~=1 && i<n^2-n)
        A(i,i+n-1)=1;
        C(i,i+n-1)=sqrt(2)*max(edges(indi,indj),edges(indi+1,indj-1));
    end
    if (mod(i,n)~=0 && i<n^2-n)
        A(i,i+n+1)=1;
        C(i,i+n+1)=sqrt(2)*max(edges(indi,indj),edges(indi+1,indj+1));
    end
end


A=A+transpose(A);
C=C+transpose(C);
center=(n^2+1)/2;
[costs,paths] = dijkstra(A,C,center);

template=zeros(n,n,n^2);

NHOOD=[0 1 0; 1 1 1; 0 1 0];
se = strel('arbitrary', NHOOD);


for i=1:n^2
    p=paths{1,i};
    % do not need to count the starting point and ending point
    for j=2:length(p)-1
        x=ceil(p(j)/n);
        y=mod(p(j), n);
        if (y==0)
            y=n;
        end
        template(x,y,i)=1;
    end
    % make 4 neighbor connected

    template(:,:,i) = imdilate(template(:,:,i), se);
    
    x1=ceil(p(1)/n);
    y1=mod(p(1), n);
    if (y1==0)
        y1=n;
    end
    template(x1,y1,i) = 1;
    
    x2=ceil(p(end)/n);
    y2=mod(p(end), n);
    if (y2==0)
        y2=n;
    end
    template(x2,y2,i) = 1;

end
