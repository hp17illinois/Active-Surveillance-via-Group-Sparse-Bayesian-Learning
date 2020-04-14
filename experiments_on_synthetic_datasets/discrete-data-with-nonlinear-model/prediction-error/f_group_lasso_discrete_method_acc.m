function [est_y,  gl_sel_set,w_inf_gl ] = f_group_lasso_discrete_method_acc(trai_x,trai_y,budget_k,nodeNum,test_x)
[w_inf_gl,gl_sel_set] = f_logistic_regression_lasso_discrete_acc(trai_x,trai_y,budget_k,nodeNum);

[duration,group_len] = size(test_x);
w_matrix = reshape(w_inf_gl,group_len,group_len);
w_matrix = w_matrix';
h_acc = 1./(1+exp(-test_x*w_matrix));
% sig_h = reshape(h_acc',duration*group_len,1);
est_y = f_trans2_01(h_acc);
end

