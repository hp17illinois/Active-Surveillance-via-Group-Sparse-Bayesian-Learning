function [mu_x,mu_x_concise,sigma_x_matrix_concise] = f_variance_aprox_embd(trai_x,trai_y,xi_concise,Sigma0,usedNum,group_len,m)

% construct a concise A
A_new =  eye(usedNum*m);
for i=1:usedNum
    for kk =1:m
        index_kk = (i-1)*m+kk;
        A_new(index_kk,index_kk) = 1/Sigma0{i}(1,1);
    end
end

index = xi_concise==0;
xi_concise(index) = 0.00001;
lam_xi_conciese = 0.5*(1./xi_concise).*(1./(1+exp(-1*xi_concise)) - 0.5);
temp_new = trai_x'*diag(lam_xi_conciese)* trai_x;

sigma_x_matrix_inv_concise = A_new +2*temp_new;
sigma_x_matrix_concise = inv(sigma_x_matrix_inv_concise);
% I = sparse(diag(ones(group_len,1)));
% sigma_x_matrix = kron(sigma_x_matrix_concise, I);

mu_x_concise =sigma_x_matrix_concise* trai_x'*(trai_y - 0.5);
mu_x = reshape(mu_x_concise',usedNum*group_len,1);
end
