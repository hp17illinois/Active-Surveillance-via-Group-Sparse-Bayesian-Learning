function [ w_inf_gl,block_ind ] = f_group_lasso_continuous_acc(attribut, label, group_len,budget_k)
%global RunTi;
EPSILON       = 1e-10;       % solution accurancy tolerance
[T,N] = size(attribut);
TN = T*N;
NN = N*N;
N = TN;
M = NN;
balance_lam = 0.2 * (N/group_len)/10;
group_num = M/group_len;
step_len = 0.0000001;
est_w = randn(M,1);
iter_num =1000;
for i=1:iter_num
%     tic;
    deno = zeros(group_num,1);
    for i_las=1:group_num
        for j_las=1:group_len
            index = (i_las-1)*group_len + j_las;
            deno(i_las) = deno(i_las) + est_w(index)^2;
        end
    end
    deno = 1./sqrt(deno);
    group_pw = balance_lam*ones(M,1);
    for i_las=1:group_num
        for j_las=1:group_len
            index = (i_las-1)*group_len + j_las;
            group_pw(index)=2*group_pw(index)*deno(i_las);
        end
    end
    A = diag(group_pw);
    
    est_w_trans = reshape(est_w,group_len,group_len);
    est_w_trans = est_w_trans';

    term_1 =  attribut' * attribut * est_w_trans;
    term_2 = attribut' * label;
    term = term_1 - term_2;
    term = term';
    term = reshape(term,group_len*group_len,1);

    
    pd = term + A*est_w;
    %pd = Phi'*(Phi*est_w-Y)+A*est_w;
    est_w_old=est_w;
    est_w = est_w - step_len*pd;
    
%     deno_t = zeros(group_num,1);
%     for i_las=1:group_num
%         for j_las=1:group_len
%             index = (i_las-1)*group_len + j_las;
%             deno_t(i_las) = deno_t(i_las) + est_w(index)^2;
%         end
%     end
%     deno_t = sqrt(deno_t);
    %recon_err = (Phi*est_w-Y)'*(Phi*est_w-Y)+ balance_lam*sum(deno_t);
    dmu = max(max(abs(est_w_old - est_w)));
    if (dmu < EPSILON)
        break;
    end;
%     toc;
    %RunTi{3} = [RunTi{3},toc];
end
group_value = zeros(group_num,1);
for i_las=1:group_num
    for j_las=1:group_len
        index = (i_las-1)*group_len + j_las;
        group_value(i_las) = group_value(i_las) + abs(est_w(index));
    end
end

[~,g_index] = sort(group_value,'descend');

est_w_tmp = zeros(group_num*group_len,1);
for node=1:budget_k
    i_las = g_index(node);
    for j_las=1:group_len
        index = (i_las-1)*group_len + j_las;
        est_w_tmp(index) = est_w(index);
    end
end

block_ind = g_index(1:budget_k);
block_ind = sort(block_ind);
block_ind = block_ind';
w_inf_gl = est_w_tmp;
end

