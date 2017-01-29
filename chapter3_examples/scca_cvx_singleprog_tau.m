function [w, e, correlation, optval, beta, mu, gamma] = scca_cvx_singleprog_tau(X, K, seed_index, sk)

% Note:
% In order to run this function, make sure that CVX toolbox is installed
% (unpack the CVX files somewhere) and run setup (call 'cvx_setup' in MATLAB
% in the CVX folder). Use the SeDuMi solver (call 'cvx_solver sedumi' after
% you have installed and run the setup). The default solver SDPT3 seems to
% run into numerical problems.

% Original description by David R. Hardoon: 
% Sparse Canonical Correlation Analysis - SCCA, is a primal-dual solver for
% the CCA problem. Given primal data of a view and a dual representation of
% the second view will provide a sparse primal weight vector (for the primal
% data) and sparse feature projection (for the dual [kernel] data)
%
% Input:  X             - Primal data of view one    [m x l]
%         K             - dual data of view two      [l x l]
%         seed_index    - Starting point for e       [1 x 1]
%         sk            - scaling factor for mu and gamma
%
% Output: w             - sparse weight vector      [1 x m]
%         e             - sparse projct vectors     [1 x l]
%         beta          - 2'nd view dual parameters [1 x l]
%         cor           - correlation value         [1 x 1]
%         optval        - Optimisation solution/s   [1 x 1]


% © 08/06/2014 Kristian Nybo, Aalto University
% kryfti@gmail.com
% 
% Original Version:
% David R. Hardoon 25/06/2007 
% http://homepage.mac.com/davidrh/
% D.Hardoon@cs.ucl.ac.uk
%
% This code is for academic purposes only.
% Commercial use is not allowed.

primal_dim = size(X,1);
N_samples = size(X,2);
tau = 0.5;

%This is how mu and gamma are set in David's SCCA2.m
Ij = eye(size(K,2));
Ij(seed_index,seed_index) = 0;
c = X*K(:,seed_index);
KK = K'*K;
d1 = 2*tau*(1-tau)*c; 
mu = sk*mean(abs(d1));
gamma = mean(abs(2*(1-tau)^2*Ij*KK(:,seed_index)));
beta = 1;

cvx_begin quiet
  variable w(primal_dim)
  variable e(N_samples)
  dual variable beta
  minimize(square_pos(norm(tau * X'*w - (1-tau)*K*e)) + mu*norm(w,1) + gamma*norm(e,1))
  beta : e >= 0
  e(seed_index) == 1
  norm(e,Inf) <= 1
cvx_end

optval = cvx_optval;
correlation = corr(X'*w, K*e);
