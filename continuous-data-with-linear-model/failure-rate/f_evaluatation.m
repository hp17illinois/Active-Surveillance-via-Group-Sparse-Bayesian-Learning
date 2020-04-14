function [ rmse_dy_pre_error, failure_rate] = f_evaluatation(est_y,test_y,hub_list,sel_set)
% dynamics prediction error, RMSE
[n,m] = size(test_y);
rmse_dy_pre_error = norm(est_y - test_y,'fro')/(n*m);
% sensor identify failure rate
len_real = length(hub_list);
C_set = intersect(hub_list, sel_set) ;
failure_rate = (len_real - length(C_set))/len_real;
end

