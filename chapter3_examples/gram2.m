function Kx = gram2(X1, X2, kernel_type, kernel_param)
% Compute the gram matrix between the sets X1 and X2
% INPUTS:
% X1:           matrix of size n1 * p
% X2:           matrix of size n2 * p
% kernel_type:  string indicating the kernel type ('linear', 'polynomial' or 'gaussian')
% kernel_param: parameter of the kernel. In the case of the polynomial
%               kernel, it corresponds to a vector containing two values
% OUTPUTS:
% Kx:           Gram matrix of size n1 * n2

%--------------------------------------------------------------------------
% © Celine Brouard, Aalto University
% celine.brouard@aalto.fi
%
% This code is for academic purposes only.
% Commercial use is not allowed.
%--------------------------------------------------------------------------

    n1 = size(X1,2);
    n2 = size(X2,2);
    
    switch kernel_type
        case 'linear'
            Kx = X1'*X2;
            
            % normalization
            Kx1 = X1'*X1;
            Kx2 = X2'*X2;
            Kx = Kx./ sqrt(repmat(diag(Kx1),1,n2) .* repmat(diag(Kx2)',n1,1));
            
        case 'polynomial'
            Kx = (X1'*X2 + kernel_param(1)).^kernel_param(2);
            
            % normalization
            Kx1 = (X1'*X1 + kernel_param(1)).^kernel_param(2);
            Kx2 = (X2'*X2 + kernel_param(1)).^kernel_param(2);
            Kx = Kx./ sqrt(repmat(diag(Kx1),1,n2) .* repmat(diag(Kx2)',n1,1));
            
        case 'gaussian'
            D = diag(X1'*X1)*ones(1,n2) + ones(n1,1)*diag(X2'*X2)' - 2*X1'*X2;
            Kx = exp(- kernel_param * D);
        otherwise
            warning('Unexpected kernel type.')
    end
    
end

