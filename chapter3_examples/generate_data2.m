function [X,Y] = generate_data2(n)
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

%n = 60;
p = 70;
q = 10;

for i = 1:p
    X(:,i) = normrnd(0,1,[n,1]);
end

y1 = X(:,3) + normrnd(0,0.01,[n 1]);
y2 = X(:,1) + normrnd(0,0.03,[n 1]);
y3 = -X(:,4) + normrnd(0,0.02,[n 1]); 

for j = 1:q-3
    Y2(:,j) = normrnd(0,1,[n,1]);
end
Y = [y1 y2 y3 Y2];

end
