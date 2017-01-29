function [za,zb,wa,wb,cc] = cca_generalised_eigenvalue_problem(X_a,X_b)

% This function solves canonical correlation analysis through 
% the generalised eigenvalue problem.

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
r = min([rank(X_a) rank(X_b)]); % dimensionality of CCA, number of relations

C_ab = (1/(n-1)) * X_a' * X_b; % empirical covariance matrix
C_ba = C_ab'; % transpose of the empirical covariance matrix
C_aa = (1/(n-1)) * X_a' * X_a; % empirical variance matrix of view a
C_bb = (1/(n-1)) * X_b' * X_b; % empirical variance matrix of view b

A = [zeros(size(C_ab,1),size(C_ba,2)) C_ab;
    C_ba zeros(size(C_ba,1),size(C_ab,2))];
B = [C_aa zeros(size(C_aa,1),size(C_bb,2));
    zeros(size(C_bb,1),size(C_aa,2)) C_bb];

[V,D] = eig(A,B); % solve the generalised eigenvalue problem
rhos = diag(D); % canonical correlations
[cc, ind] = sort(rhos,'descend'); % sort in descending order
vecs = V(:,ind); 
wa = vecs(1:size(X_a,2),1:r);
wb = vecs(size(X_a,2)+1:end,1:r);

za = normc(X_a * wa); % on the surface of the unit ball
zb = normc(X_b * wb); % on the surface of the unit ball



end

