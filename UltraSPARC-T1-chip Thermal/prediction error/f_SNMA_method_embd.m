function [ est_y,bsbl_sel_set,new_networks,gamma_value ] = f_SNMA_method_embd( attribut_embd,label,k,test_x,m)
[sampleTimeLength,nodeNum]=size(attribut_embd);
m=m+1;
nodeNum = nodeNum/m;
   
Result = f_SNMA_continuous_concise_embd(attribut_embd,label,nodeNum,sampleTimeLength,k,m);

est_W = Result.x;
network = reshape(est_W,nodeNum,nodeNum*m );
new_networks = network';
est_y = test_x * new_networks;
%new_networks = network * sum(sum(gt_networks))/(sum(sum(network)));
bsbl_sel_set = Result.gamma_used';
gamma_value = Result.gamma_est';
end

