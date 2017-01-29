function [za,zb,wa,wb,cc,T,U,S,V] = cca_svd(X_a,X_b)

% This function solves canonical correlation analysis using 
% the singular value decomposition.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       X_a: the data matrix of the 1st view (view a)
%       X_b: the data matrix of the 2nd view (view b)

% Output
%       za: the image of the position wa (a.k.a. canonical variate)
%       zb: the image of the position wb (a.k.a. canonical variate)
%       wa: the position in the data space of view a
%       wb: the position in the data space of view b
%       cc: the canonical correlation (cosine of the enclosing angle)
%--------------------------------------------------------------------------
% © 19/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

n = size(X_a,1); % number of training observations
C_ab = (1/(n-1)) * X_a' * X_b; % empirical covariance matrix
C_aa = (1/(n-1)) * X_a' * X_a; % empirical variance matrix of view a
C_bb = (1/(n-1)) * X_b' * X_b; % empirical variance matrix of view b

T = C_aa^(-1/2) * C_ab * C_bb^(-1/2);
[U,S,V] = svd(T);
sigma = diag(S);
[cc, inds] = sort(sigma,'descend');
U = U(:,inds);
V = V(:,inds);

wa = C_aa^(-1/2) * U;
wb = C_bb^(-1/2) * V;

za = normc(X_a * wa);
zb = normc(X_b * wb);


end

