function [ rmse_dy_pre_error] = f_evaluatation(est_y,test_y)
% dynamics prediction error, RMSE
[n,m] = size(test_y);
rmse_dy_pre_error = norm(est_y - test_y,'fro')/(n*m);
end

