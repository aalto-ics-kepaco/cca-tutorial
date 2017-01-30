function [X_a,X_b] = generate_data(n)
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

% normally distributed variables in the view a
a1 = normrnd(0,1,[n,1]);
a2 = normrnd(0,1,[n,1]);
a3 = normrnd(0,1,[n 1]);
a4 = normrnd(0,1,[n 1]);

% linearly related variables in the view b plus normal noise
b1 = a3 + normrnd(0,0.2,[n 1]);
b2 = a1 + normrnd(0,0.4,[n 1]);
b3 = -a4 + normrnd(0,0.3,[n 1]); 

X_a = [a1 a2 a3 a4];
X_b = [b1 b2 b3];

end
