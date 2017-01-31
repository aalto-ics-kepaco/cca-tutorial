function [final, c1_opt, c2_opt] = crossv_kernel(X,Y,F)

% This function performs repeated cross-validation to determine the optimal
% regularisation parameters of kernel canonical correlation analysis.

% Uurtio et al. A Tutorial on Canonical Correlation Methods. 2017.
%--------------------------------------------------------------------------
% Input 
%       X: the data matrix of the 1st view
%       Y: the data matrix of the 2nd view
%       F: number of folds
% Output
%       final: the test canonical correlations at all c1 and c2 values
%       c1_opt: the optimal c1
%       c2_opt: the optimal c2
%--------------------------------------------------------------------------
% © 30/01/2017 Viivi Uurtio, Aalto University
% viivi.uurtio@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

pats = 10; % number of relations, patterns
c1 = 0.2:0.1:1.5; % values to be tested for c1
c2 = 0.2:0.1:1.5; % values to be tested for c1

for rep = 1:20 % repetitions
    [~,indices] = stratified_cv(size(X,1), F);
    for j = 1:size(c1,2)
        for k = 1:size(c2,2)
            for i = 1:F
                
                test = indices == i;
                train = indices ~= i;
                A = zscore(X(train,:));
                B = zscore(Y(train,:));
                n = size(A,1);
                
                [Ka,siga] = gaussK(A,'median');
                [Kb,sigb] = gaussK(B,'median');
                Ka = centering(Ka); Kb = centering(Kb);
                
                Kar = Ka + c1(j) * eye(size(Ka));
                Kbr = Kb + c2(k) * eye(size(Kb));
                
                O = [zeros(n) Ka*Kb;
                    Kb*Ka zeros(n)];
                P = [Kar^2 zeros(n);
                    zeros(n) Kbr^2];
                
                [Vg,Dg] = eig(O,P);
                rhosg = diag(Dg);
                [~, ind] = sort(rhosg,'descend');
                vecsg = Vg(:,ind);
                
                wag = vecsg(1:n,1:pats);
                wbg = vecsg(n+1:end,1:pats);
                
                %sag = normc(Ka * wag);
                %sbg = normc(Kb * wbg);
                
                Atest = zscore(X(test,:));
                Btest = zscore(Y(test,:));
                
                Katest = gram2(Atest', A', 'gaussian', siga);
                Kbtest = gram2(Btest', B', 'gaussian', sigb);
                
                %Katest = centering(Katest); Kb = centering(Kbtest);
                
                za = normc(Katest * wag(:,1));
                zb = normc(Kbtest * wbg(:,1));
                
                test_corr(i,1) = za'*zb;
                
            end
            mean_corr(j,k,rep) = mean(test_corr);
        end
    end
    final = mean(mean_corr,3);
end


[max_val, ind] = max(final(:));
[r, c] = ind2sub(size(final),ind);
mean_opt = [r c];
c1_opt = c1(r); c2_opt = c2(c);
[c1_opt c2_opt]

end

