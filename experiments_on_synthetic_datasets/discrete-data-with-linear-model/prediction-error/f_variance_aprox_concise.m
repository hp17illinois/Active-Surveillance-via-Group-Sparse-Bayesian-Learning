function [mu_x,mu_x_concise,sigma_x_matrix_concise_set] = f_variance_aprox_concise(trai_x,trai_y,xi,Sigma0,nodeNum,usedNum,sampleTimeLength)
% construct a concise A
A_new =  eye(usedNum);
for i=1:usedNum
    A_new(i,i) = 1/Sigma0{i}(1,1);
end

% index = xi_concise==0;
% xi_concise(index) = 0.00001;
% lam_xi_conciese = 0.5*(1./xi_concise).*(1./(1+exp(-1*xi_concise)) - 0.5);

index = xi ==0;
xi(index) = 0.00001;
lam_xi = 0.5*(1./xi).*(1./(1+exp(-1*xi)) - 0.5);

sigma_x_matrix_concise_set = cell(nodeNum,1);
mu_x_tmp = zeros(usedNum,nodeNum);
for i = 1: nodeNum
    index = i:nodeNum:nodeNum*sampleTimeLength;
    lam_xi_tmp = lam_xi(index);
    temp_new = trai_x'*diag(lam_xi_tmp)* trai_x;
    sigma_x_matrix_inv_concise = A_new +2*temp_new;
    sigma_x_matrix_concise = inv(sigma_x_matrix_inv_concise);
    sigma_x_matrix_concise_set{i} = sigma_x_matrix_concise;
    mu_x_tmp(:,i) =sigma_x_matrix_concise* trai_x'*(trai_y(:,i) - 0.5);
end
    
% temp_new = trai_x'*diag(lam_xi_conciese)* trai_x;
% 
% sigma_x_matrix_inv_concise = A_new +2*temp_new;
% 
% sigma_x_matrix_concise = inv(sigma_x_matrix_inv_concise);
% % I = sparse(diag(ones(group_len,1)));
% % sigma_x_matrix = kron(sigma_x_matrix_concise, I);
% % 
% % mu_x_concise =sigma_x_matrix_concise* trai_x'*(trai_y - 0.5);
% % % mu_x_concise = pinv(pinv(trai_y - 0.5) *  pinv(trai_x')* sigma_x_matrix_inv_concise); 
% mu_x = reshape(mu_x_concise',usedNum*nodeNum,1);
mu_x_concise = mu_x_tmp;
mu_x = reshape(mu_x_tmp',usedNum*nodeNum,1);
end
