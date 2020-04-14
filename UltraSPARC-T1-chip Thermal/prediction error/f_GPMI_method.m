function [ est_y,MI_sel_set ] = f_GPMI_method( trai_x,trai_y,budget_k,test_x )
%%% training
[~, nodeNum] = size(trai_x);
full_data = horzcat(trai_x,trai_y);
Sigma = cov(full_data);

v_set = [1:nodeNum]; % observation
u_set = [nodeNum+1:2*nodeNum]; % prediction
a_set = []; % selected set
c_set = v_set; % candidature

%%% submodual greed strategy 
tic
for i = 1:budget_k
    [a_set,c_set] = f_GPMI_submodual(Sigma,u_set,a_set,c_set);
end
MI_sel_set = sort(a_set); % sensor locations
x = test_x(:,MI_sel_set);
mean_fun_x = mean(trai_x(:,MI_sel_set));
mean_fun_y = mean(trai_y);
%%% GP dynamics prediction
est_y = f_GPMI_dynamicsPrediction(x,Sigma,MI_sel_set,u_set,mean_fun_x,mean_fun_y);
end

