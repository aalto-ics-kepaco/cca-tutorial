function [X,Y] = generate_scca_data(n)

% This function generates normally distributed data containing the 
% specific relations studied in the examples.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       n: number of observations

% Output
%       X_a: The data matrix of view a
%       X_b: The data matrix of view a
%--------------------------------------------------------------------------
% © 30/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

%n = 50;
p = 100;
q = 150;

for i = 1:p
    X(:,i) = normrnd(0,1,[n,1]);
end

for j = 1:q
    Y(:,j) = normrnd(0,1,[n,1]);
end




end
