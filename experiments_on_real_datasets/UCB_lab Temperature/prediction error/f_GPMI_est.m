function [ est_y_mean ] = f_GPMI_est( est_node_id, t, mean_fun_x,mean_fun_y, x,Sig_aa,Sig_ya  )

mean_y = mean_fun_y(est_node_id);
obser_x = x(t,:)';

est_y_mean = mean_y + Sig_ya(est_node_id,:) *pinv(Sig_aa) * (obser_x - mean_fun_x');

end

