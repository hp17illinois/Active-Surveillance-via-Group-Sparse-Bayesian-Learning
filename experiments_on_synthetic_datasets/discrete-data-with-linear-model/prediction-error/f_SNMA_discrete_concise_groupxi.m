function Result = f_SNMA_discrete_concise_groupxi(trai_x_0,trai_y,nodeNum,sampleTimeLength,sen_k)
% Default Parameter Values for Any Cases
EPSILON       = 1e-20;       % solution accurancy tolerance
MAX_ITERS     = 500;        % maximum iterations
PRINT         = 0;          % don't show progress information
lambda = 1e-3;
PRUNE_GAMMA = 1e-300;

%% Initialization
[~,M] = size(trai_x_0);
M = M*M;
blkLen = nodeNum;
blkStartLoc = [1:blkLen:M];
blkStartLoc0 = blkStartLoc;

blkLenList = blkLen * ones(1,nodeNum);
maxLen = nodeNum;
for k = 1 : nodeNum
    Sigma0{k} = eye(blkLenList(k));
end

gamma = ones(nodeNum,1);
keep_list = [1:nodeNum]';
usedNum = length(keep_list);
mu_x = zeros(M,1);
xi_concise = 100*ones(sampleTimeLength,1);
count = 0;
init_flag = 1;

trai_x=trai_x_0;
%% Iteration
while (length(gamma)>sen_k)
    if init_flag == 1
        init_flag = 0;
    else
        [~,I] = min(gamma);
        gamma(I)=0;
        count = 0;
        %=========== Prune weights as their \gamma go to zero ==============
        index = find(gamma > 1e-6);
        usedNum = length(index);
        keep_list = keep_list(index);
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        % prune gamma and associated nodes in Sigma0 and data
        gamma = gamma(index);
        temp = Sigma0;
        
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        temp = [];
        for k = 1 : usedNum
            temp = [temp, trai_x_0(:,keep_list(k))];
        end
        trai_x = temp;
    end
    
    while (1)
        if count > MAX_ITERS
            break;
        end
        count = count +1;
        %=================== Compute new weights =================
        mu_old = mu_x;
        Sigma_x = [];
        Cov_x = [];
      
        [mu_x,mu_x_concise,sigma_x_matrix_concise] = f_variance_aprox_group_xi(trai_x,trai_y,xi_concise, Sigma0,usedNum,maxLen);

        xi_concise = [];
        for i=1:sampleTimeLength
            temp_1 = trai_x(i,:)*sigma_x_matrix_concise*trai_x(i,:)';
            temp_total = 0;
            xi_tmp = zeros(1,nodeNum);
            for j = 1: nodeNum
                temp_3 = trai_x(i,:)*mu_x_concise(:,j);
                temp_3 = temp_3*temp_3;
                temp_2 = sqrt(temp_1+temp_3);
                xi_tmp(j) = temp_2;
                temp_total = temp_total + temp_2;
            end
			xi_concise = [xi_concise;temp_total/nodeNum]; % average xi
        end
        
        for i = 1 : usedNum
            Sigma_x{i}  = eye(maxLen)*sigma_x_matrix_concise(i,i);
            Cov_x{i} = Sigma_x{i} + mu_x_concise(i,:)' * mu_x_concise(i,:);
        end
        currentLoc = 0;
        for i =  1 : usedNum
            gamma(i) = trace( Cov_x{i})/size(Cov_x{i},1);
            [tmp_n,~]=size(Sigma0{i});
            tmp = ones(tmp_n,1)*gamma(i);
            Sigma0{i} = diag(tmp);
        end

        % ================= Check stopping conditions, eyc. ==============
        if (size(mu_x) == size(mu_old))
            dmu = max(max(abs(mu_old - mu_x)));
            if (dmu < EPSILON)  break;  end;
        end;
    end;
end


gamma_used = sort(keep_list);
gamma_est = zeros(nodeNum,1);
gamma_est(keep_list,1) = gamma;

x = zeros(M,1);
currentLoc = 0;
for i = 1 : usedNum
    currentLen = size(Sigma0{i},1);
    currentLoc = currentLoc + 1;
    seg = currentLoc : 1 : currentLoc + currentLen - 1;
    realLocs = blkStartLoc0(keep_list(i)) : blkStartLoc0(keep_list(i))+currentLen-1;
    x( realLocs ) = mu_x( seg );
    currentLoc = seg(end);
end

Result.x = x;
Result.gamma_used = gamma_used;
Result.count = count;
return;

