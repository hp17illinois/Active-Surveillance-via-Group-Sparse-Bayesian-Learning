function [ MNEP_est_y, MNEP_sel_set ] = f_MNEP_PCA_method( trai_x,trai_y,budget_k,test_x )
[~, nodeNum] = size(trai_x);
data_X = horzcat(trai_x,trai_y);
[Phi_frse] = pca(data_X);
data_X = data_X';
[~,k] = size(Phi_frse); 
k=budget_k; % the max k
Phi_frse_k = Phi_frse(:,1:k); % all k eigen vectors

% laten_repr = Phi_frse_k' * data_X; %reconstruction
% data_recons = Phi_frse_k * laten_repr;
% 
% Phi_frse_k_s = Phi_frse_k(1:budget_k,:); %random select and reconstruction
% data_recons_s_p2 =  Phi_frse_k * pinv(transpose(Phi_frse_k_s)*Phi_frse_k_s)*transpose(Phi_frse_k_s) * data_X(1:budget_k,:);

forb_list = [nodeNum+1:2*nodeNum]; % prediction, cannot deploy sensor;
feasible_list = [1:nodeNum]; % the locations wherer sensors can be deployed
SensorPosition = f_MNEP_greed_forb(Phi_frse_k,budget_k,feasible_list);
Phi_frse_k_s = zeros(budget_k, k);
data_S=[];
test_x=test_x';
for i = 1:budget_k
    Phi_frse_k_s(i,:) = Phi_frse_k(SensorPosition(i),:);
    data_S(i,:) = test_x(SensorPosition(i),:);
end
data_prediction_s_p_fra =  Phi_frse_k *  pinv(transpose(Phi_frse_k_s)*Phi_frse_k_s) * transpose(Phi_frse_k_s)  * data_S;
MNEP_sel_set = sort(SensorPosition');
MNEP_est_y = data_prediction_s_p_fra(forb_list,:)';
end





