function [ est_y, gl_sel_set,w_inf_gl ] = f_group_lasso_discrete_method( trai_Phi,trai_Y,budget_k,nodeNum,test_Phi )
[w_inf_gl,gl_sel_set] = f_logistic_regression_lasso(trai_Y,trai_Phi,budget_k,nodeNum);
h = -1 * test_Phi * w_inf_gl;
sig_h = 1./(1+exp(h));
est_y = f_trans2_01(sig_h);
end

