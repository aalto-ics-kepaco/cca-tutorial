function L = bartlett_lawley_statistic(n,k,p,q,r)

% This function computes the Bartlett-Lawley statistic for Bartlett's
% sequential test procedure.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       n: number of observations, sample size
%       k: number of canonical correlations to be tested
%       p: number of variables in the 1st view
%       q: number of variables in the 2nd view
%       r: a vector containing the canonical correlations

% Output
%       L: the Bartlett-Lawley statistic to be compared against the
%       chi-squared distribution with (p-k)(q-k) degrees of freedom.
%--------------------------------------------------------------------------
% © 19/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

psum = 1;
for j = k+1:min(p,q)
    psum = psum * (1-r(j)^2);
end

L = -(n - k - 1/2 * (p + q + 1) + sum(r.^(-2))) * log(psum);

end

