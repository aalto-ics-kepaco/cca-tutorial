%% Numerical Examples of Chapter 2

% This m file contains scripts to reproduce the results of the numerical
% examples of the chapter 2 of the paper "A Tutorial on Canonical 
% Correlation Methods".

%--------------------------------------------------------------------------
% © 29/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

%% Example 2.1. Solving CCA Through the Standard Eigenvalue Problem
clear
rng(0) % initialise the random number generator to reproduce the results

n = 60; % number of observations
[X_a,X_b] = generate_data(n); % observations from a random normal multivariate distribution

[~,indices] = stratified_cv(size(X_a,1), 3); % partition into training and test sets
test = indices == 1; % all 1s constitute the test set (1/3)
train = indices ~= test; % all 2s and 3s constitute the training set (2/3)

X_atest = X_a(test,:); X_btest = X_b(test,:); % select the test observations, needed for Example 2.5
X_a = X_a(train,:); X_b = X_b(train,:); % select the training observations 
X_a = zscore(X_a); X_b = zscore(X_b); % standardise the training variables

[za_sta, zb_sta, wa_sta, wb_sta, cc_sta, ev_sta] = cca_standard_eigenvalue_problem(X_a,X_b);

%% Example 2.2. Solving CCA Through the Generalised Eigenvalue Problem

[za_gen, zb_gen, wa_gen, wb_gen, cc_gen] = cca_generalised_eigenvalue_problem(X_a,X_b);

%% Example 2.3. Solving CCA Using the SVD

[za_svd, zb_svd, wa_svd, wb_svd, cc_svd] = cca_svd(X_a, X_b);

%% Example 2.4. Bartlett's sequential test procedure

n = size(X_a,1); % sample size
p = size(X_a,2); % number of variables in view a
q = size(X_b,2); % number of variables in view b

r = cc_sta; % the canonical correlations
k = 0:1:min(p,q)-1; % number of canonical correlations
for i = 1:size(k,2)
    L(i) = bartlett_lawley_statistic(n,k(i),p,q,r); % the statistics
    prob(i) = 1 - chi2cdf(L(i),(p-k(i))*(q-k(i)));
    x(i) = chi2inv(0.99,(p-k(i))*(q-k(i))); % critical values from the chi-squared distribution
end

%% Example 2.5. The generalisability of the canonical correlation model

za_test = normc(zscore(X_atest) * wa_sta);
zb_test = normc(zscore(X_btest) * wb_sta);

% test canonical correlations
[za_test(:,1)'*zb_test(:,1); 
    za_test(:,2)'*zb_test(:,2); 
    za_test(:,3)'*zb_test(:,3)]



