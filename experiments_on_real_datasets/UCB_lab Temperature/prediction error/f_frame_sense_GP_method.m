function [ FsGP_est_y, FsGP_sel_set ] = f_frame_sense_GP_method( trai_x,trai_y,budget_k,test_x )
[sample_duration, nodeNum] = size(test_x);
full_data = horzcat(trai_x,trai_y);
Sigma = cov(full_data);
s_set = [1:nodeNum]; % observation location, where we can place sensors
u_set = [nodeNum+1:2*nodeNum]; % prediction locations
Sigma_ss = Sigma(s_set,s_set);
Sigma_us = Sigma(u_set,s_set);
A = Sigma_ss*pinv(Sigma_us);
SensorPosition = f_frame_sense_greed_forb(A,budget_k,[]);
FsGP_sel_set = SensorPosition;

A_s = A(SensorPosition,:);
data_x_s = test_x(:,SensorPosition); 

mean_x = mean(trai_x);
mean_x = mean_x(SensorPosition);
mean_matrix_x=repmat(mean_x,sample_duration,1);
mean_y = mean(trai_y);
mean_matrix_y=repmat(mean_y,sample_duration,1);


est_y = pinv(transpose(A_s)*A_s)*transpose(A_s)*transpose(data_x_s-mean_matrix_x)+transpose(mean_matrix_y);
est_y = transpose(est_y);
FsGP_est_y =est_y;

end





