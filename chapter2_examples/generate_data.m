function [X_a,X_b] = generate_data(n)

% This function generates normally distributed data containing the 
% specific relations studied in the examples.

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
