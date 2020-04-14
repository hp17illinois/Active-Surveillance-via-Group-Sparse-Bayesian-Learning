function [ H ] = f_GPMI_Gaussian_entropy( Sigma )
[k,~] = size(Sigma);
 H = 0.5 * log( (2*pi*exp(1))^k * det(Sigma)); % at lesat O(N^3) complexity
 
% if det(Sigma) <= 0
%     H =   -1.0000e+73;
% else
%     H = 0.5 * log( (2*pi*exp(1))^k * det(Sigma)); % at lesat O(N^3) complexity
% end
    


end

