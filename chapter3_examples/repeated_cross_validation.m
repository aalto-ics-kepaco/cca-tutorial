function [c1_opt,c2_opt,final] = repeated_cross_validation(X_a,X_b,c1,c2,reps)
% This function performs repeated cross-validation to determine the optimal
% regularisation parameters.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       X_a: the data matrix of the 1st view
%       X_b: the data matrix of the 2nd view
%       c1: a vector containing the values to be tested
%       c2: a vector containing the values to be tested
%       reps: number of repetitions of the cross-validation
% Output
%       c1_opt: the optimal c1
%       c2_opt: the optimal c2
%       final: the test canonical correlations at all c1 and c2 values
%--------------------------------------------------------------------------
% � 30/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

F = 5;
for rep = 1:reps
    [~,indices] = stratified_cv(size(X_a,1), F);
    for j = 1:size(c1,2)
        for k = 1:size(c2,2)            
            for i = 1:5            
                test = indices == i;
                train = indices ~= i;
                
                % standardise
                trainX_a = zscore(X_a(train,:));
                trainX_b = zscore(X_b(train,:));
                
                % covariance matrix
                n = size(trainX_a,1);
                C_ab = (1/n) * trainX_a' * trainX_b;
                C_ba = C_ab';
                C_aa = (1/n) * trainX_a' * trainX_a;
                C_bb = (1/n) * trainX_b' * trainX_b;
                
                % regularise
                C_aa = C_aa + c1(j) * eye(size(C_aa));
                C_bb = C_bb + c2(k) * eye(size(C_bb));
                
                M = inv(C_bb) * C_ba * inv(C_aa) * C_ab;
                [eigvectors,eigvalues] = eig(M) ;
                rhos = diag(sqrt(eigvalues));
                
                % sort the eigenvalues and eigenvectors
                [eig_sorted, ind] = sort(rhos,'descend');
                wb = eigvectors(:,ind);
		wb = wb(:,1); % we only consider the first eigenvector
                
                wa = zeros(size(trainX_a,2),1);
                for i = 1:size(wb,2) 
                    wa(:,i) = (inv(C_aa) * C_ab * wb(:,i)) / eig_sorted(i);
                end                
                za = normc(zscore(X_a(test,:)) * wa(:,1));
                zb = normc(zscore(X_b(test,:)) * wb(:,1));
                
                test_corr(i,1) = za'*zb;
            end        
            mean_corr(j,k,rep) = mean(test_corr);  
        end
    end   
end
final = mean(mean_corr,3);

[~, ind] = max(final(:));
[r, c] = ind2sub(size(final),ind);
mean_opt = [r, c];
c1_opt = c1(r); c2_opt = c2(c);
[c1_opt c2_opt]


end

