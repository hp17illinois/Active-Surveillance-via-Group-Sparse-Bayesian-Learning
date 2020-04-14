function [ est_y, gl_sel_set,w_inf_gl ] = f_group_lasso_method(Phi,Y,budget_k,test_set_x,nodeNum )

tic;
[ w_inf_gl,block_ind ] = f_group_lasso_continuous( Phi,Y, budget_k,nodeNum);
toc;

network = reshape(w_inf_gl,nodeNum,nodeNum );
network = network';
est_y = test_set_x * network;
gl_sel_set = block_ind;

end

