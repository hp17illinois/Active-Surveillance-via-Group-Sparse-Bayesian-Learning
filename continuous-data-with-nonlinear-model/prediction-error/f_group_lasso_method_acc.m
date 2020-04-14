function [ est_y, gl_sel_set,w_inf_gl] = f_group_lasso_method_acc(attribut,label,budget_k,test_set_x )

[~,nodeNum]=size(attribut);

tic;
[ w_inf_gl,block_ind ] = f_group_lasso_continuous_acc(attribut, label, nodeNum,budget_k);
toc;

network = reshape(w_inf_gl,nodeNum,nodeNum );
network = network';
est_y = test_set_x * network;
gl_sel_set = block_ind;

end

