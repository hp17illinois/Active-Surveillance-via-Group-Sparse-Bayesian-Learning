clear all;
warning off all
%% step1: env setting %%
tran_test_splite_point = 160;
tran_test_splite_point = 80;
budget_k_list = [5,10,15,20,25,30,35,40,45]; % sensor number according to budget.
[~, len_k] = size(budget_k_list);

%% step2:load epidemic dynamics %%
[attributs,attributsEmbd , label, duration, nodeNum] = loadRealData_lab_embd();

trai_x = attributs(1:tran_test_splite_point,:);
trai_x_embd = attributsEmbd(1:tran_test_splite_point,:);
trai_y = label(1:tran_test_splite_point,:);

test_x = attributs(tran_test_splite_point+1:end,:);
test_x_embd = attributsEmbd(tran_test_splite_point+1:end,:);
test_y = label(tran_test_splite_point+1:end,:);
%% step3:algorithm comparison %%

rmse_matrix = [];
for i = 1:len_k
    rmse_tmp = [];
    budget_k = budget_k_list(i); % sensor number according to budget.
    %=== Gaussian processing MI ===%
    tic
    [MI_est_y,MI_sel_set ] = f_GPMI_method( trai_x,trai_y,budget_k,test_x );
    [ rmse_mi] = f_evaluatation(MI_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_mi];
    % mse(1,cv_i) = mse_mi;
    toc
    %=== SNMA ===%
    tic
    [ SNMA_est_y, SNMA_sel_set,est_network,gamma_value] = f_SNMA_method( trai_x,trai_y,budget_k,test_x );
    [ rmse_SNMA] = f_evaluatation(SNMA_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_SNMA];
    toc
    %=== SNMA_embd ===%
    tic
    NumBasicFunc =1;
    [ SNMA_embd_est_y, SNMA_embd_sel_set,~] = f_SNMA_method_embd( trai_x_embd,trai_y,budget_k,test_x_embd,NumBasicFunc);
    [ rmse_SNMA] = f_evaluatation(SNMA_embd_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_SNMA];
    toc
    %===  Group_lasso ===%
    tic
    [ GL_est_y, Gl_sel_set,w_inf_gl ] = f_group_lasso_method_acc( trai_x,trai_y,budget_k,test_x );
    [rmse_gl] = f_evaluatation(GL_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_gl];
    toc
    %===  Frame sense based on PCA ===%
    tic
    [ Fs_PCA_est_y, Fs_PCA_sel_set ] = f_frame_sense_PCA_method( trai_x,trai_y,budget_k,test_x );
    [rmse_Fs_PCA] = f_evaluatation(Fs_PCA_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_Fs_PCA];
    toc
    %===  MNEP based on PCA ===%
    tic
    [ MNEP_PCA_est_y, MNEP_PCA_sel_set ] = f_MNEP_PCA_method( trai_x,trai_y,budget_k,test_x );
    [rmse_MNEP_PCA] = f_evaluatation(MNEP_PCA_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_MNEP_PCA];
    toc
    %===  MPME based on PCA ===%
    tic
    [ MPME_PCA_est_y, MPME_PCA_sel_set ] = f_MPME_PCA_method(trai_x,trai_y,budget_k,test_x );
    [rmse_MPME_PCA] = f_evaluatation(MPME_PCA_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_MPME_PCA];
    toc
    %===  Frame sense based on GP ===%
    tic
    [ Fs_GP_est_y, Fs_GP_sel_set ] = f_frame_sense_GP_method( trai_x,trai_y,budget_k,test_x );
    [rmse_Fs_GP] = f_evaluatation(Fs_GP_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_Fs_GP];
    toc
    %===  MNEP based on GP ===%
    tic
    [ MNEP_GP_est_y, MNEP_GP_sel_set ] = f_MNEP_GP_method( trai_x,trai_y,budget_k,test_x );
    [rmse_MNEP_GP] = f_evaluatation(MNEP_GP_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_MNEP_GP];
    toc
    %===  MPME based on GP ===%
    tic
    [ MPME_GP_est_y, MPME_GP_sel_set ] = f_MPME_GP_method(trai_x,trai_y,budget_k,test_x );
    [rmse_MPME_GP] = f_evaluatation(MPME_GP_est_y,test_y);
    rmse_tmp = [rmse_tmp,rmse_MPME_GP];
    toc
    rmse_matrix = [rmse_matrix;rmse_tmp];
end


