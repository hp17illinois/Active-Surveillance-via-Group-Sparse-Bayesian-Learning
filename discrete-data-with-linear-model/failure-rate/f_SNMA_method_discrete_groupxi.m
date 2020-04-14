function [ est_y,bsbl_sel_set,new_networks ] = f_SNMA_method_discrete_groupxi(trai_x,trai_y,sen_k,test_x)

[sampleTimeLength,nodeNum] = size(trai_x);
% Result = f_SNMA_discrete_concise_working(trai_Phi,trai_Y,trai_x,trai_y,nodeNum, sampleTimeLength,sen_k);
Result = f_SNMA_discrete_concise_groupxi(trai_x,trai_y,nodeNum,sampleTimeLength,sen_k);
% Result = BSBL_EM_discre(trai_Phi,trai_Y,trai_x,trai_y,nodeNum,sen_k);

est_W = Result.x;
network = reshape(est_W,nodeNum,nodeNum );
new_networks = network';
%new_networks = network * sum(sum(gt_networks))/(sum(sum(network)));


h_acc = 1./(1+exp(-test_x*new_networks));
% sig_h = reshape(h_acc',duration*group_len,1);
est_y = f_trans2_01(h_acc);

bsbl_sel_set = Result.gamma_used';
end

