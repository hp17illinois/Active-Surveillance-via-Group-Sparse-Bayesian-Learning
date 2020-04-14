function updat_w = f_for_up_pd_acc(inf_w,trai_x,trai_y,A)
[~,group_len] = size(trai_x);
w_matrix = reshape(inf_w,group_len,group_len);
w_matrix = w_matrix';

mu_acc = 1./(1+exp(-trai_x*w_matrix));

%mu = 1./(1+exp(-x*inf_w));
% pd_w = x' * (y-mu);

pd_w_acc=trai_x'*(trai_y - mu_acc);
pd_w = reshape(pd_w_acc',group_len*group_len,1);
pd_w = pd_w - A*inf_w;

% V = mu.*(1-mu);
% V = diag(V);
% H = -x'*V*x-A;
% updat_w = H\pd_w;

updat_w = -pd_w;