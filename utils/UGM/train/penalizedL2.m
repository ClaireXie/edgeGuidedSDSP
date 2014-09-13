function [nll,g,H] = penalizedL2(w,gradFunc,lambda,varargin)
% [nll,g,H] = penalizedL2(w,gradFunc,lambda,varargin)
% Adds L2-penalization to a loss function
% (you can use this instead of always adding it to the loss function code)

if nargout <= 1
    [nll] = gradFunc(w,varargin{:});
elseif nargout == 2
    [nll,g] = gradFunc(w,varargin{:});
else
    [nll,g,H] = gradFunc(w,varargin{:});
end

nll = nll+sum(lambda.*(w.^2));

if nargout > 1
    g = g + 2*lambda.*w;
end

if nargout > 2
    if isscalar(lambda)
        H = H + 2*lambda*eye(length(w));
    else
        H = H + diag(2*lambda);
    end
end