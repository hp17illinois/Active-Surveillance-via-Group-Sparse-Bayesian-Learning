function [w_inf_gl,block_ind] = f_logistic_regression_lasso_discrete_acc(trai_x,trai_y,budget_k,group_len)
[n,m]=size(trai_x);
banl_beta = 1 * (n/m);
m = m*m;
step_len = 0.01;
group_num = m/group_len;
EPSILON       = 1e-5;       % solution accurancy tolerance
inf_w = rand(m,1);
iter_num =2000;
J_pre = 10000000000;
count_error = 0;
thresh_error = 10;% 
for i =1:iter_num
    deno = zeros(group_num,1);
    for i_las=1:group_num
        for j_las=1:group_len
            index = (i_las-1)*group_len + j_las;
            deno(i_las) = deno(i_las) + inf_w(index)^2;
        end
    end
    deno = 1./sqrt(deno);  %group lasso panety
    group_pw = banl_beta*ones(m,1);
    for i_las=1:group_num
        for j_las=1:group_len
            index = (i_las-1)*group_len + j_las;
            group_pw(index)=2*group_pw(index)*deno(i_las);
        end
    end
    A = diag(group_pw);
    updat_w = f_for_up_pd_acc(inf_w,trai_x,trai_y,A);
    inf_w_old=inf_w;
    inf_w = inf_w - step_len * updat_w;
    dmu = max(max(abs(inf_w_old - inf_w)));
    if (dmu < EPSILON)
        break;
    end;
    J = f_group_lasso_discrete_costFunction_acc( inf_w,trai_x,trai_y,A );
    if J > J_pre
        count_error = count_error +1;
    else
        count_error = 0;
    end
    if count_error > thresh_error
        break;
    end
    J_pre = J;
end

group_value = zeros(group_num,1);
for i_las=1:group_num
    for j_las=1:group_len
        index = (i_las-1)*group_len + j_las;
        group_value(i_las) = group_value(i_las) + abs(inf_w(index));
    end
end
[~,g_index] = sort(group_value,'descend');
est_w_tmp = zeros(group_num*group_len,1);
for node=1:budget_k
    i_las = g_index(node);
    for j_las=1:group_len
        index = (i_las-1)*group_len + j_las;
        est_w_tmp(index) = inf_w(index);
    end
end
block_ind = g_index(1:budget_k);
block_ind = sort(block_ind);
block_ind = block_ind';
w_inf_gl = est_w_tmp;
end
