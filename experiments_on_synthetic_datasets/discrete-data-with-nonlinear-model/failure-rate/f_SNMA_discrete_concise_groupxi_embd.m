function Result = f_SNMA_discrete_concise_groupxi_embd(trai_x_0,trai_y,nodeNum,sampleTimeLength,sen_k,m)
% Default Parameter Values for Any Cases
EPSILON       = 1e-10;       % solution accurancy tolerance
MAX_ITERS     = 400;        % maximum iterations
PRINT         = 0;          % don't show progress information
lambda = 1e-3;
PRUNE_GAMMA = 1e-300;

%% Initialization
M = nodeNum*nodeNum*m;
blkLen = nodeNum*m;
blkStartLoc = [1:blkLen:M];
blkStartLoc0 = blkStartLoc;
blkLenList = blkLen * ones(1,nodeNum);

maxLen = nodeNum*m;

Sigma0_new = eye(nodeNum*m);
for k = 1 :nodeNum
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
        index = find(gamma > 0);
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
        
        %construct new attributes
        temp = [];
        for k = 1 : usedNum
            ind_tem = keep_list(k);
            real_index_start = (ind_tem-1)*m+1;
            real_index_end = m*ind_tem;
            for real_index = real_index_start : real_index_end
                temp = [temp, trai_x_0(:,real_index)];
            end
        end
        trai_x = temp;
    end
    
    while (1)
        count = count + 1;
        %=================== Compute new weights =================
        mu_old = mu_x;
        Sigma_x = [];
        Cov_x = [];
        
        [mu_x,mu_x_concise,sigma_x_matrix_concise] = f_variance_aprox_embd(trai_x,trai_y,xi_concise,Sigma0,usedNum,maxLen,m);

        xi_concise = [];
        for i=1:sampleTimeLength
            temp_1 = trai_x(i,:)*sigma_x_matrix_concise*trai_x(i,:)';
            temp_total = 0;
            for j = 1: nodeNum
                temp_2 = trai_x(i,:)*mu_x_concise(:,j)*mu_x_concise(:,j)'*trai_x(i,:)';
                temp_total = temp_total + sqrt(temp_1+temp_2);
            end
            xi_concise = [xi_concise;temp_total/nodeNum]; %average xi
        end

        for i = 1 : usedNum
            S_tmp = [];
            mu_x_tmp = [];
            for kk = 1:m
                real_ind = (i-1)*m+kk;
                tmp = ones(nodeNum,1)*sigma_x_matrix_concise(real_ind,real_ind);
                S_tmp = [S_tmp;tmp];
                mu_x_tmp = [mu_x_tmp,mu_x_concise(real_ind,:)];
            end
            Sigma_x{i}  = diag(S_tmp);
            Cov_x{i} = Sigma_x{i} + mu_x_tmp * mu_x_tmp';
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
        if count > MAX_ITERS
            break;
        end
        
    end;
end


gamma_used = sort(keep_list);

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

