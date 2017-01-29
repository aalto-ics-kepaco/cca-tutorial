function [Kc] = centering(K)
% center the kernel matrices
ell = size(K,1);
D = sum(K) / ell;
E = sum(D) / ell;
J = ones(ell,1) * D;
Kc = K - J - J' + E * ones(ell,ell);


end

