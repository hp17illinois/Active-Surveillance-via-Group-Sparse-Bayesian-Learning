function J = f_group_lasso_discrete_costFunction_acc( inf_w,trai_x,trai_y,A )
[duration,group_len] = size(trai_x);
w_matrix = reshape(inf_w,group_len,group_len);
w_matrix = w_matrix';
h_acc = 1./(1+exp(-trai_x*w_matrix));
h = reshape(h_acc',duration*group_len,1);

y = reshape(trai_y',duration*group_len,1);
J = (h-y)'*(h-y) + inf_w'*A*inf_w;

end


