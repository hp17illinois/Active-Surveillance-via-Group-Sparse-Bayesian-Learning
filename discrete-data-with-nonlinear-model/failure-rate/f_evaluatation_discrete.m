function [ precision, failure_rate] = f_evaluatation_discrete(est_y,test_y,hub_list,sel_set)
% dynamics prediction error, RMSE

 [n,m] = size(test_y);
% test_y = reshape(test_y',n*m,1);
diff = est_y == test_y;
precision = sum(sum(diff))/(n*m);

% sensor identify failure rate
len_real = length(hub_list);
C_set = intersect(hub_list, sel_set) ;
failure_rate = (len_real - length(C_set))/len_real;
end

