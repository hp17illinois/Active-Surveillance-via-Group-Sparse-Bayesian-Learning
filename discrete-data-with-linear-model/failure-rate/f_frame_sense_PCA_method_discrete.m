function [ FsPCA_est_y, FsPCA_sel_set ] = f_frame_sense_PCA_method_discrete( trai_x,trai_y,budget_k,test_x )
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
SensorPosition = f_frame_sense_greed_modified(Phi_frse_k,budget_k,forb_list);
Phi_frse_k_s = zeros(budget_k, k);
data_S=[];
test_x=test_x';
for i = 1:budget_k
    Phi_frse_k_s(i,:) = Phi_frse_k(SensorPosition(i),:);
    data_S(i,:) = test_x(SensorPosition(i),:);
end
data_prediction_s_p_fra =  Phi_frse_k *  pinv(transpose(Phi_frse_k_s)*Phi_frse_k_s) * transpose(Phi_frse_k_s)  * data_S;
FsPCA_sel_set = SensorPosition;
FsPCA_est_y = data_prediction_s_p_fra(forb_list,:)';

ind_1 = FsPCA_est_y > 0.5;
ind_0 = FsPCA_est_y <= 0.5;
FsPCA_est_y(ind_1) = 1;
FsPCA_est_y(ind_0) = 0;

end





