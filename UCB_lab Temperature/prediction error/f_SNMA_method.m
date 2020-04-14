function [ est_y,bsbl_sel_set,new_networks,gamma_value] = f_SNMA_method( attribut,label,k,test_x )
[sampleTimeLength,nodeNum]=size(attribut);
     
Result = f_SNMA_continuous_concise(attribut,label,nodeNum,sampleTimeLength,k);

est_W = Result.x;
network = reshape(est_W,nodeNum,nodeNum );
new_networks = network';
est_y = test_x * new_networks;
%new_networks = network * sum(sum(gt_networks))/(sum(sum(network)));
bsbl_sel_set = Result.gamma_used';
gamma_value = Result.gamma_est';
end

