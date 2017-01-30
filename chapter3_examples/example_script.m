%% Examples in Chapter 3

% This m file contains scripts to reproduce the results of the numerical
% examples of the chapter 3 of the paper "A Tutorial on Canonical 
% Correlation Methods".

%--------------------------------------------------------------------------
% © 30/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------
%% Example 3.1. Regularised CCA
clear % clear workspace
rng(0) % initialise the random number generator to reproduce the results

[X_a,X_b] = generate_data2(60);  % observations from a random normal multivariate distribution

c1 = 0.01:0.01:1; % testing range for regularisation parameter c1
c2 = 0; % testing range for regularisation parameter c2
reps = 50; % number of repetitions for cross-validation
[c1_opt,c2_opt,final] = repeated_cross_validation(X_a,X_b,c1,c2,reps); % find optimal regularisation parameters

X_a = zscore(X_a); X_b = zscore(X_b); % standardise
[za_sep,zb_sep,wa_sep,wb_sep,cc_sep,ev_sep] = cca_standard_regularised(X_a,X_b,c1_opt,c2_opt); % solve CCA through standard eigenvalue problem
[za_gep,zb_gep,wa_gep,wb_gep,cc_gep,ev_gep] = cca_generalised_regularised(X_a,X_b,c1_opt,c2_opt); % or through the generalised eigenvalue problem

%% Example 3.2. Kernel CCA
clear % clear workspace
rng(0) % initialise the random number generator to reproduce the results

[X_a,X_b] = generate_kernel_data2(150);  % observations from a random normal multivariate distribution
X_a = zscore(X_a); X_b = zscore(X_b); % standardise

[Ka,siga] = gaussK(X_a,'median'); % apply the median trick to determine the width and kernelise
[Kb,sigb] = gaussK(X_b,'median'); % apply the median trick to determine the width and kernelise

Ka = centering(Ka); Kb = centering(Kb); % center the kernels
rels = 3; % number of relations of interest

[mean_corr, c1_opt, c2_opt] = crossv_kernel(X_a,X_b,5); % find the optimal regularisation parameters

[za_sep,zb_sep,alpha_sep,beta_sep,cc_sep,ev_sep] = cca_standard_kernel(Ka,Kb,c1_opt,c2_opt,rels); % solve by sep
[za_gep,zb_gep,alpha_gep,beta_gep,cc_gep,ev_gep] = cca_generalised_kernel(Ka,Kb,c1_opt,c2_opt,rels); % solve by gep

[cc_sep' cc_gep'] % take a look at the canonical correlations

% which relations are extracted by which pairs of images 
pat1 = [abs(corr(exp(X_a(:,3)),za_gep))' abs(corr(X_b(:,1),zb_gep))']
pat2 = [abs(corr(X_a(:,1).^3,za_gep))' abs(corr(X_b(:,2),zb_gep))']
pat3 = [abs(corr(-X_a(:,4),za_gep))' abs(corr(X_b(:,3),zb_gep))']

%% Example 3.3. Kernel CCA using PGSO

clear
rng(0) % initialise the random number generator to reproduce the results

[X_a,X_b] = generate_kernel_data2(10000); % generate large data matrices
X_a = zscore(X_a); X_b = zscore(X_b); % standardise

[Ka,siga] = gaussK(X_a,'median'); % apply the median trick to determine the width and kernelise
[Kb,sigb] = gaussK(X_b,'median'); % apply the median trick to determine the width and kernelise

Ka = centering(Ka); Kb = centering(Kb); % center the kernels
rels = 3; % number of relations of interest

[nalpha, nbeta, r, Kx, Ky] = kcanonca_reg_ver2(Ka,Kb,0.1,0.5,0,2);

zah = normc(Ka * nalpha); 
zbh = normc(Kb * nbeta); 

zah = zah(:,[1:rels]);
zbh = zbh(:,[1:rels]);

for i = 1:3
    corr_genh(i) = zah(:,i)'*zbh(:,i);
end
[cans, inds] = sort(corr_genh,'descend');
zah = zah(:,inds); zbh = zbh(:,inds);

pat1 = [abs(corr(exp(X_a(:,3)),zah))' abs(corr(X_b(:,1),zbh))'];
pat2 = [abs(corr(X_a(:,1).^3,zah))' abs(corr(X_b(:,2),zbh))'];
pat3 = [abs(corr(-X_a(:,4),zah))' abs(corr(X_b(:,3),zbh))'];


%% Example 3.4. Sparse CCA Through PMD
clear
load pmd_result.mat % load the positions wa, wb obtained by pmd in R
%load sparse_data.mat
rng(0)
[X_a,X_b] = generate_data_pmd(50); % the data matrices

za_pmd = normc(zscore(X_a) * wa); % on the surface of the unit ball
zb_pmd = normc(zscore(X_b) * wb); % on the surface of the unit ball

% the first three canonical correlations
[za_pmd(:,1)'*zb_pmd(:,1);
za_pmd(:,2)'*zb_pmd(:,2);
za_pmd(:,3)'*zb_pmd(:,3)]

%% Example 3.5. Primal-Dual Sparse CCA
% To run this example, the cvx toolbox is required.
clear
rng(0)

[X_a,X_b] = generate_scca_data(50); % generate data
X_a = zscore(X_a); X_b = zscore(X_b); % standardise

X = X_a;
[K,sig] = gaussK(X_b,'median'); % Gaussian kernel, width by median trick
[K] = centering(K); % center the kernel

% checking which k results in the minimum objective value
k = size(X_a,1);
for i = 1:k
    [wa, e, corval, optval, beta, mu, gamma] = scca_cvx_singleprog_tau(X',K,i,1);
    was(:,i) = wa;
    es(:,i) = e;
    corvals(i) = corval;
    optvals(i) = optval;
    mus(i) = mu;
    gammas(i) = gamma;
end

% put very small values to zero
was(was < 1e-6) = 0;
es(es < 1e-6) = 0;

[val, ind] = min(optvals);
min_opt = optvals(ind); 
optcor = corvals(ind); % the optimal canonical correlation
W = was(:,ind); % variables in view a
E = es(:,ind); % observations in view b




