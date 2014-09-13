function filtered = bilateralOMA(D)

Dmi = min(D(:));
Dn = D-Dmi;
Dma = max(Dn(:));
Dn = Dn/Dma;

%filtered = bfilter2(Dn,5,[1.5 0.01]); % less smoothing
filtered = bfilter2(Dn,5,[1.5 0.05]);
%filtered = bfilter2(Dn,5,[3 0.1]); % more smoothing 0.1

filtered = (filtered*Dma)+Dmi;