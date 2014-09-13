
%Useful function to perform computations in log domain.
% M - matrix
% i = 1 or 2
% returns log(sum(exp(M),i)), but carefully computed in log domain to
% avoid numerical issues
function res=logsumExp(M,i)
    scl =   max(M,[],i);
    M   =   bsxfun(@minus,M,scl);
    res =   log(sum(exp(M),i))+scl;
    res =   res(:);