function [za,zb,wa,wb,cc,ev] = cca_generalised_regularised(X_a,X_b,c1,c2)

% This function solves regularised canonical correlation analysis through 
% the generalised eigenvalue problem.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       X_a: the data matrix of the 1st view (view a)
%       X_b: the data matrix of the 2nd view (view b)
%       c1: the regularisation parameter for the 1st view
%       c2: the regularisation parameter for the 2nd view
% Output
%       za: the image of the position wa (a.k.a. canonical variate)
%       zb: the image of the position wb (a.k.a. canonical variate)
%       wa: the position in the data space of view a
%       wb: the position in the data space of view b
%       cc: the canonical correlation (cosine of the enclosing angle)
%       ev: the generalised eigenvalues
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

% regularise
C_aa_reg = C_aa + c1 * eye(size(C_aa));
C_bb_reg = C_bb + c2 * eye(size(C_bb));

% the matrix pencil
A = [zeros(size(C_ab,1),size(C_ba,2)) C_ab;
    C_ba zeros(size(C_ba,1),size(C_ab,2))];
B = [C_aa_reg zeros(size(C_aa,1),size(C_bb,2));
    zeros(size(C_bb,1),size(C_aa,2)) C_bb_reg];

[V,D] = eig(A,B); % solve the generalised eigenvalue problem
rhos = diag(D); % canonical correlations
[ev, ind] = sort(rhos,'descend'); % sort in descending order
vecs = V(:,ind); 
wa = vecs(1:size(X_a,2),1:r);
wb = vecs(size(X_a,2)+1:end,1:r);

za = normc(X_a * wa); % on the surface of the unit ball
zb = normc(X_b * wb); % on the surface of the unit ball

for i = 1:r
    cc(i) = za(:,i)'*zb(:,i);
end
[cc, ix] = sort(cc,'descend');
za = za(:,ix); zb = zb(:,ix);



end

