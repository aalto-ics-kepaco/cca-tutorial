function [za,zb,alpha,beta,cc,ev] = cca_standard_kernel(Ka,Kb,c1,c2,rels)

% This function solves kernel canonical correlation analysis through 
% the standard eigenvalue problem.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       Ka: the Gram matrix of the 1st view (view a)
%       Kb: the Gram matrix of the 2nd view (view b)
%       c1: the regularisation parameter for the 1st view
%       c2: the regularisation parameter for the 2nd view
%       rels: number of relations of interest

% Output
%       za: the image of the position wa (a.k.a. canonical variate)
%       zb: the image of the position wb (a.k.a. canonical variate)
%       alpha: the position in the data space of view a
%       beta: the position in the data space of view b
%       cc: the canonical correlation (cosine of the enclosing angle)
%       ev: the square roots of the eigenvalues
%--------------------------------------------------------------------------
% © 19/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

% regularise
Kar = (Ka + c1 * eye(size(Ka)))^2;
Kbr = (Kb + c2 * eye(size(Kb)))^2;

M = inv(Kbr) * Kb * Ka* inv(Kar) * Ka * Kb; % the canonical correlation matrix

[eigvectors,eigvalues] = eig(M) ; % solve the standard eigenvalue problem
rhos = diag(sqrt(eigvalues)); % canonical correlations
[ev, ind] = sort(rhos,'descend'); % sort in descending order
wb = eigvectors(:,ind); % the positions in the data space of view b
beta = wb(:,1:rels);

alpha = zeros(size(Ka,1),rels);
for i = 1:size(beta,2)
    alpha(:,i) = (inv(Kar) * Ka * Kb * beta(:,i)) / ev(i); 
end
 
za = normc(Ka * alpha);
zb = normc(Kb * beta);

for i = 1:rels
    cc(i) = za(:,i)'*zb(:,i);
end
[cc, ix] = sort(cc,'descend');
za = za(:,ix); zb = zb(:,ix);


end