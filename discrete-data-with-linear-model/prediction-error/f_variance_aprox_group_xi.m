% function [mu_x,mu_x_concise,sigma_x_matrix_concise] = f_variance_aprox(trai_x,trai_y,xi_concise,xi,Sigma0,usedNum,group_len)
function [mu_x,mu_x_concise,sigma_x_matrix_concise] = f_variance_aprox_group_xi(trai_x,trai_y,xi_concise,Sigma0,usedNum,group_len)
% construct a concise A
A_new =  eye(usedNum);
for i=1:usedNum
    A_new(i,i) = 1/Sigma0{i}(1,1);
end

index = xi_concise==0;
xi_concise(index) = 0.00001;
lam_xi_conciese = 0.5*(1./xi_concise).*(1./(1+exp(-1*xi_concise)) - 0.5);

% index = xi ==0;
% xi(index) = 0.00001;
% lam_xi = 0.5*(1./xi).*(1./(1+exp(-1*xi)) - 0.5);

temp_new = trai_x'*diag(lam_xi_conciese)* trai_x;

sigma_x_matrix_inv_concise = A_new +2*temp_new;

sigma_x_matrix_concise = inv(sigma_x_matrix_inv_concise);
% I = sparse(diag(ones(group_len,1)));
% sigma_x_matrix = kron(sigma_x_matrix_concise, I);

mu_x_concise =sigma_x_matrix_concise* trai_x'*(trai_y - 0.5);
% mu_x_concise = pinv(pinv(trai_y - 0.5) *  pinv(trai_x')* sigma_x_matrix_inv_concise); 

mu_x = reshape(mu_x_concise',usedNum*group_len,1);
end
