function [cv,indices] = stratified_cv(y, nfold)
% This function partitions the data into training and test sets.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       y: vector of labels
%       nfold: number of folds
% Output
%       cv: struct from cvpartition
%       indices: the indices of the partitioning
%--------------------------------------------------------------------------
% © 19/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

cv = cvpartition(y,'kfold',nfold);
indices = zeros(size(y));
for q = 1:nfold
    indices(cv.test(q)) = q;
end





end