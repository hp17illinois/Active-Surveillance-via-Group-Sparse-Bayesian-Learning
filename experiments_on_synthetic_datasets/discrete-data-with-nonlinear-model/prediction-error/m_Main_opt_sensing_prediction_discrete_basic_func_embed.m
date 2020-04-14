
clear all;

budget_k =20; % sentinel/sensor number according to budget.
%% step1:generate a network %%
k_hub = 20;    % pre-given number of important nodes
nodeNum=50;   % number of nodes
if (budget_k > nodeNum) or  (k_hub > nodeNum)
    error('too many hubs or sensors');
end
sparsity = 0.1;
NumBasicFunc=1; % the number of ground truth basic function
[network, Hub_list] = f_netGen_given_k_disc_basic_fun_embd(nodeNum, k_hub, sparsity,NumBasicFunc); % scaled network and sorted hub list

%% step2:generate epidemic dynamics %%
trai_sampleTimeLength = 70;  % sample number in training dateset
test_sampleTimeLength = 70;   % sample number in testing dateset
BER = 0.15; % noise level

periodInterval = 5;
[trai_x,trai_x_embd, trai_y , Wgen] = f_dynamicsGen_compac_discretet_poly_fun_embd(trai_sampleTimeLength,network,NumBasicFunc,periodInterval,BER);
[test_x,test_x_embd, test_y , ~]    = f_dynamicsGen_compac_discretet_poly_fun_embd(test_sampleTimeLength,network,NumBasicFunc,periodInterval,BER);
%% step3:algorithm comparison %%
budget_k_list = [2,4,6,8,10,12,14,16,18,20]; % sensor number according to budget.
[~, len_k] = size(budget_k_list);
rmse_matrix = [];
for i = 1:len_k
    rmse_tmp = [];
    budget_k = budget_k_list(i); % sensor number according to budget.
    % %=== Gaussian processing MI ===%
    tic
    [MI_est_y,MI_sel_set ] = f_GPMI_method_discrete( trai_x,trai_y,budget_k,test_x );
    [ precision_MI, failure_rate_MI] = f_evaluatation_discrete(MI_est_y,test_y,Hub_list,MI_sel_set);
    rmse_tmp = [rmse_tmp,precision_MI];
    toc
    %=== SNMA ===%
    tic
    [ SNMA_est_y_groupxi, SNMA_sel_set_groupxi,est_network] = f_SNMA_method_discrete_groupxi_embd( trai_x_embd,trai_y,budget_k,test_x_embd,NumBasicFunc);
    [ precision_SNMA_groupxi_embd, failure_rate_SNMA_groupxi_embd] = f_evaluatation_discrete(SNMA_est_y_groupxi,test_y,Hub_list,SNMA_sel_set_groupxi);
    toc
    %===  Group_lasso classification ===%
    tic
    [ GL_est_y, Gl_sel_set,w_inf_gl ] = f_group_lasso_discrete_method_acc(trai_x,trai_y,budget_k,nodeNum,test_x );
    [ precision_GL, failure_rate_GL] = f_evaluatation_discrete(GL_est_y,test_y,Hub_list,Gl_sel_set);
    rmse_tmp = [rmse_tmp,precision_GL];
    toc
    %===  Frame sense based on PCA ===%
    tic
    [ Fs_PCA_est_y, Fs_PCA_sel_set ] = f_frame_sense_PCA_method_discrete( trai_x,trai_y,budget_k,test_x );
    [precision_Fs_PCA, failure_rate_Fs_PCA] = f_evaluatation_discrete(Fs_PCA_est_y,test_y,Hub_list,Fs_PCA_sel_set);
    rmse_tmp = [rmse_tmp,precision_Fs_PCA];
    toc
    %===  MNEP based on PCA ===%
    tic
    [ MNEP_PCA_est_y, MNEP_PCA_sel_set ] = f_MNEP_PCA_method_discrete( trai_x,trai_y,budget_k,test_x );
    [precision_MNEP_PCA, failure_rate_MNEP_PCA] = f_evaluatation_discrete(MNEP_PCA_est_y,test_y,Hub_list,MNEP_PCA_sel_set);
    rmse_tmp = [rmse_tmp,precision_MNEP_PCA];
    toc
    %===  MPME based on PCA ===%
    tic
    [ MPME_PCA_est_y, MPME_PCA_sel_set ] = f_MPME_PCA_method_discrete(trai_x,trai_y,budget_k,test_x );
    [precision_MPME_PCA, failure_rate_MPME_PCA] = f_evaluatation_discrete(MPME_PCA_est_y,test_y,Hub_list,MPME_PCA_sel_set);
    rmse_tmp = [rmse_tmp,precision_MPME_PCA];
    toc
    %===  Frame sense based on GP ===%
    tic
    [ Fs_GP_est_y, Fs_GP_sel_set ] = f_frame_sense_GP_method_discrete( trai_x,trai_y,budget_k,test_x );
    [precision_Fs_GP, failure_rate_Fs_GP] = f_evaluatation_discrete(Fs_GP_est_y,test_y,Hub_list,Fs_GP_sel_set);
    rmse_tmp = [rmse_tmp,precision_Fs_GP];
    toc
    %===  MNEP based on GP ===%
    tic
    [ MNEP_GP_est_y, MNEP_GP_sel_set ] = f_MNEP_GP_method_discrete( trai_x,trai_y,budget_k,test_x );
    [precision_MNEP_GP, failure_rate_MNEP_GP] = f_evaluatation_discrete(MNEP_GP_est_y,test_y,Hub_list,MNEP_GP_sel_set);
    rmse_tmp = [rmse_tmp,precision_MNEP_GP];
    toc
    %===  MPME based on GP ===%
    tic
    [ MPME_GP_est_y, MPME_GP_sel_set ] = f_MPME_GP_method_discrete(trai_x,trai_y,budget_k,test_x );
    [precision_MPME_GP, failure_rate_MPME_GP] = f_evaluatation_discrete(MPME_GP_est_y,test_y,Hub_list,MPME_GP_sel_set);
    rmse_tmp = [rmse_tmp,precision_MPME_GP];
    toc
    rmse_matrix = [rmse_matrix;rmse_tmp];
end
avg_prediction_error = mean(rmse_matrix);
fprintf('GP-MI avg_prediction_error: %f \n', avg_prediction_error(1));
fprintf('SNMA avg_prediction_error: %f \n',  avg_prediction_error(2));
fprintf('Group_lasso avg_prediction_error: %f \n',  avg_prediction_error(3));
fprintf('Fs_PCA avg_prediction_error: %f \n',  avg_prediction_error(4));
fprintf('MNEP_PCA avg_prediction_error: %f \n',  avg_prediction_error(5));
fprintf('MPME_PCA avg_prediction_error: %f \n',  avg_prediction_error(6));
fprintf('Fs_GP avg_prediction_error: %f \n',  avg_prediction_error(7));
fprintf('MNEP_GP avg_prediction_error: %f \n',  avg_prediction_error(8));
fprintf('MPME_GP avg_prediction_error: %f \n',  avg_prediction_error(9));



