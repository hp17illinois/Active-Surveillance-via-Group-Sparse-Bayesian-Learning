function est_y = f_GPMI_dynamicsPrediction(x,Sigma,Sel_set,u_set,mean_fun_x,mean_fun_y)

est_node_number = length(u_set);
[est_duration,~] = size(x);
est_y = zeros(est_duration,est_node_number);


Sig_aa = Sigma(Sel_set,Sel_set);
Sig_ya = Sigma(u_set:end, Sel_set);

for t = 1: est_duration
    for est_node_id = 1:est_node_number
        [ est_y_mean ] = f_GPMI_est( est_node_id, t, mean_fun_x,mean_fun_y , x,Sig_aa,Sig_ya);
        est_y(t,est_node_id) = est_y_mean;
    end
end

end

