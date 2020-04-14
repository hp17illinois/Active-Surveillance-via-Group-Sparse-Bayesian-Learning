function [mu_x,sigma_x_matrix] = logistic_regression(y,Phi,xi,Sigma0,group_len)

[n,m]=size(Phi);
k=length(Sigma0);
A = eye(k*group_len); %total sigma_0 inversed
for i=1:k
    A( (i-1)*group_len+1:i*group_len,(i-1)*group_len+1:i*group_len ) = inv(Sigma0{i});
end

% temp = ones(n,1)*0.5;
% sig_xi = 1./(1+exp(-xi));
% Lamda = -0.5 * (sig_xi-temp)./xi;
% Lamda = diag(Lamda);
% sigma_x_matrix = A + 2*Phi'*Lamda*Phi;
% sigma_x_matrix = inv(sigma_x_matrix);

temp = zeros(m,m);
lam_xi_list = zeros(m,1);
for i =1:n
    if xi(i) ~= 0
       lam_xi = 0.5*(1/xi(i))*(1/(1+exp(-xi(i))) - 0.5);
    else
        lam_xi = 0.5;
    end
    lam_xi_list(i) = lam_xi;
    Phi_list{i} = lam_xi * Phi(i,:)'* Phi(i,:);
    temp= temp+lam_xi* Phi(i,:)'* Phi(i,:);
end

sigma_x_matrix_inv = A + 2*temp;
sigma_x_matrix = inv(sigma_x_matrix_inv);

temp = zeros(m,1);
for i = 1:n
   temp = temp + (y(i)-0.5)*Phi(i,:)';
end
mu_x = sigma_x_matrix * temp;

end
