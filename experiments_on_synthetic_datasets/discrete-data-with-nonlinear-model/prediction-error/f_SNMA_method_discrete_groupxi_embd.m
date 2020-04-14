function [ est_y,bsbl_sel_set,new_networks ] = f_SNMA_method_discrete_groupxi_embd(trai_x_embd,trai_y,sen_k,test_x_embd,m)

[sampleTimeLength,nodeNum] = size(trai_x_embd);
m=m+1;
nodeNum = nodeNum/m;
Result = f_SNMA_discrete_concise_groupxi_embd(trai_x_embd,trai_y,nodeNum,sampleTimeLength,sen_k,m);

est_W = Result.x;
network = reshape(est_W,nodeNum,nodeNum*m );
new_networks = network';

h_acc = 1./(1+exp(-test_x_embd*new_networks));
% sig_h = reshape(h_acc',duration*group_len,1);
est_y = f_trans2_01(h_acc);

bsbl_sel_set = Result.gamma_used';
end

