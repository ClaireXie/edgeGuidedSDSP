function C=constructGraphWeight(row, col, edges, n, C0)

C = zeros(n^2);

for i = 1:length(row)
    indi1 = ceil(row(i)/n);
    indj1 = mod(row(i)-1,n)+1;
    
    indi2 = ceil(col(i)/n);
    indj2 = mod(col(i)-1,n)+1;
    
    C(row(i),col(i)) = C0(row(i),col(i))*max(edges(indi1,indj1),edges(indi2,indj2));
end
