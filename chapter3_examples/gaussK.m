function [K,sig] = gaussK(X,type,sigma)
%%return kernel K and bandwidth found by median trick 

switch type
    case 'median'
        K = distanceMatrix(X);
        sig = median(K(:));
        K = exp(- (K.^2)./(2*sig.^2));
    case 'none'
        K = distanceMatrix(X);
        K = exp(- (K.^2)./(2*sigma.^2));
        sig = sigma;
end

end 

function D = distanceMatrix(X)

N = size(X,1);

XX = sum(X.*X,2);
XX1 = repmat(XX,1,N);
XX2 = repmat(XX',N,1);

D = XX1+XX2-2*(X*X');
D(D<0) = 0;
D = sqrt(D);
end