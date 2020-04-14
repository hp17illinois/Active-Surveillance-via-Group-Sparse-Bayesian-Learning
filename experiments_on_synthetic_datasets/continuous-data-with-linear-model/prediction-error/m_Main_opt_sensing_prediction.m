clear all;

%% step1:generate a network %%
k_hub = 20;    % pre-given number of important nodes
nodeNum= 50;   % number of nodes
% if (budget_k > nodeNum) or  (k_hub > nodeNum)
%     error('too many hubs or sensors');
% end
sparsity = 0.1;
[network, Hub_list] = f_netGen_given_k(nodeNum, k_hub, sparsity); % scaled network and sorted hub list

%% step2:generate epidemic dynamics %%
periodInterval = 2;  % length of a continuous epidemic period
trai_sampleTimeLength = 80;  % sample number in training dateset
test_sampleTimeLength = 80;   % sample number in testing dateset

[trai_x, trai_y , Wgen] = f_dynamicsGen_compact(trai_sampleTimeLength,network,periodInterval);
[test_x, test_y , ~]    = f_dynamicsGen_compact(test_sampleTimeLength,network,periodInterval);
SNR = 10;
[trai_y] = f_add_noise_continu(trai_y,SNR); % add continuous niose according to given SNR

%% step3:algorithm comparison %%
budget_k_list = [2,4,6,8,10,12,14,16,18,20]; % sentinel/sensor number according to budget.
[~, len_k] = size(budget_k_list);
rmse_matrix = []; %prediction error
for i = 1:len_k %len_k
    rmse_tmp = [];
    budget_k = budget_k_list(i); % sensor number according to budget.
    %=== Gaussian processing MI ===%
    tic
    [MI_est_y,MI_sel_set ] = f_GPMI_method( trai_x,trai_y,budget_k,test_x );
    [ rmse_mi, failure_mi] = f_evaluatation(MI_est_y,test_y,Hub_list,MI_sel_set);
    rmse_tmp = [rmse_tmp,rmse_mi];
    toc
    %=== SNMA ===%
    tic
    [ SNMA_est_y, SNMA_sel_set,est_network] = f_SNMA_method( trai_x,trai_y,budget_k,test_x );
    [ rmse_SNMA, failure_SNMA] = f_evaluatation(SNMA_est_y,test_y,Hub_list,SNMA_sel_set);
    rmse_tmp = [rmse_tmp,rmse_SNMA];
    toc
    %===  Group_lasso ===%
    tic
    [ GL_est_y, Gl_sel_set,w_inf_gl ] = f_group_lasso_method_acc( trai_x,trai_y,budget_k,test_x );
    [rmse_gl, failure_gl] = f_evaluatation(GL_est_y,test_y,Hub_list,Gl_sel_set);
    rmse_tmp = [rmse_tmp,rmse_gl];
    toc
    %===  Frame sense based on PCA ===%
    tic
    [ Fs_PCA_est_y, Fs_PCA_sel_set ] = f_frame_sense_PCA_method( trai_x,trai_y,budget_k,test_x );
    [rmse_Fs_PCA, failure_Fs_PCA] = f_evaluatation(Fs_PCA_est_y,test_y,Hub_list,Fs_PCA_sel_set);
    rmse_tmp = [rmse_tmp,rmse_Fs_PCA];
    toc
    %===  MNEP based on PCA ===%
    tic
    [ MNEP_PCA_est_y, MNEP_PCA_sel_set ] = f_MNEP_PCA_method( trai_x,trai_y,budget_k,test_x );
    [rmse_MNEP_PCA, failure_MNEP_PCA] = f_evaluatation(MNEP_PCA_est_y,test_y,Hub_list,MNEP_PCA_sel_set);
    rmse_tmp = [rmse_tmp,rmse_MNEP_PCA];
    toc
    %===  MPME based on PCA ===%
    tic
    [ MPME_PCA_est_y, MPME_PCA_sel_set ] = f_MPME_PCA_method(trai_x,trai_y,budget_k,test_x );
    [rmse_MPME_PCA, failure_MPME_PCA] = f_evaluatation(MPME_PCA_est_y,test_y,Hub_list,MPME_PCA_sel_set);
    rmse_tmp = [rmse_tmp,rmse_MPME_PCA];
    toc
    %===  Frame sense based on GP ===%
    tic
    [ Fs_GP_est_y, Fs_GP_sel_set ] = f_frame_sense_GP_method( trai_x,trai_y,budget_k,test_x );
    [rmse_Fs_GP, failure_Fs_GP] = f_evaluatation(Fs_GP_est_y,test_y,Hub_list,Fs_GP_sel_set);
    rmse_tmp = [rmse_tmp,rmse_Fs_GP];
    toc
    %===  MNEP based on GP ===%
    tic
    [ MNEP_GP_est_y, MNEP_GP_sel_set ] = f_MNEP_GP_method( trai_x,trai_y,budget_k,test_x );
    [rmse_MNEP_GP, failure_MNEP_GP] = f_evaluatation(MNEP_GP_est_y,test_y,Hub_list,MNEP_GP_sel_set);
    rmse_tmp = [rmse_tmp,rmse_MNEP_GP];
    toc
    %===  MPME based on GP ===%
    tic
    [ MPME_GP_est_y, MPME_GP_sel_set ] = f_MPME_GP_method(trai_x,trai_y,budget_k,test_x );
    [rmse_MPME_GP, failure_MPME_GP] = f_evaluatation(MPME_GP_est_y,test_y,Hub_list,MPME_GP_sel_set);
    rmse_tmp = [rmse_tmp,rmse_MPME_GP];
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



