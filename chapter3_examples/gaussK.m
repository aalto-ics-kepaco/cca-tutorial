function [K,sig] = gaussK(X,type,sigma)
% This function returns the Gaussian gram matrix and sets the kernel width
% using the median trick if requested.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       X: the data matrix (rows - observations, columns - variables)
%       type: 'median' for median trick, 'none' for width selected by user
%       sigma: if 'none', set the width
% Output
%       K: The gram matrix
%       sig: the width of the kernel
%--------------------------------------------------------------------------
% © 30/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

switch type
    case 'median'
        K = distanceMatrix(X);
        sig = median(K(:));
        K = exp(- (K.^2)./(2*sig.^2));
    case 'none'
        K = distanceMatrix(X);
        K = exp(- (K.^2)./(2*sigma.^2));
        sig = sigma;
end

end 

function D = distanceMatrix(X)

N = size(X,1);

XX = sum(X.*X,2);
XX1 = repmat(XX,1,N);
XX2 = repmat(XX',N,1);

D = XX1+XX2-2*(X*X');
D(D<0) = 0;
D = sqrt(D);
end