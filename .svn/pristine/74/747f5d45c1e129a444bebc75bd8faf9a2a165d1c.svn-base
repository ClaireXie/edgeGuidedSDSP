function D = depthToDisp( ip, gt, scaleFact )

D = 1./(gt*scaleFact);
minGt = min(D(:));
D = D-minGt;
maxGt = max(D(:));
D = D/maxGt;

D = (1./((ip*maxGt)+minGt))/scaleFact;

end


% convert from disp to depth
%D = 1./D;
%D = D-min(D(:));
%D = D/max(D(:));