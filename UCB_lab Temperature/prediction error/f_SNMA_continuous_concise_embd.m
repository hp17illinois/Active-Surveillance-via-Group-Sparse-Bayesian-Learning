function Result = f_SNMA_continuous_concise_embd(attribut,label, nodeNum,sampleTimeLength,sen_k,m)
%global RunTi;
% Default Parameter Values for Any Cases
EPSILON       = 1e-20;       % solution accurancy tolerance
MAX_ITERS     =1000;        % maximum iterations
lambda = 1e-3;
lambda_new = 1e-3;
gamma_epis = 1e-6;
PRINT =  0;

%% Initialization
TN = sampleTimeLength * nodeNum*m;
NN = nodeNum * nodeNum*m;

attribut0 = attribut;
blkLen = nodeNum*m;
blkStartLoc = [1:blkLen:NN];
blkStartLoc0 = blkStartLoc;
blkLenList = blkLen * ones(1,nodeNum);

maxLen = nodeNum*m;

Sigma0_new = eye(nodeNum*m);
for k = 1 : nodeNum
    Sigma0{k} = eye(blkLenList(k));
end

gamma = ones(nodeNum,1);
keep_list = [1:nodeNum]';
usedNum = length(keep_list);
mu_x = zeros(NN,1);
count = 0;
init_flag = 1;
%% Iteration
%sen_k =nodeNum-sen_k ;
while (length(gamma)>sen_k)
    if init_flag == 1
        init_flag = 0;
    else
        [~,I] = min(gamma);
        gamma(I)=0;
        count = 0;
        %=========== Prune weights as their hyperparameters go to zero ==============
        index = find(gamma > gamma_epis);
        while (length(index) < sen_k)
            gamma_epis = gamma_epis *0.5;
            index = find(gamma > gamma_epis);
        end
            
        usedNum = length(index);
        usedNum
        keep_list = keep_list(index);
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        
        % prune gamma and associated components in Sigma0
        gamma = gamma(index);
        temp = Sigma0;
        
        gamma_embd = [];
        for k = 1:usedNum
            for kk = 1:m
                gamma_embd = [gamma_embd,gamma(k)];
            end
        end
        Sigma0_new = diag(gamma_embd);
        
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new attributes
        temp = [];
        for k = 1 : usedNum
            ind_tem = keep_list(k);
            real_index_start = (ind_tem-1)*m+1;
            real_index_end = m*ind_tem;
            for real_index = real_index_start : real_index_end
                temp = [temp, attribut0(:,real_index)];
            end
        end
        attribut = temp;
    end
    while (1)
        count = count + 1;
    
        %=================== Compute new weights =================
        mu_old = mu_x;
        
        PhiBPhi_new = attribut*Sigma0_new*attribut';
        H_new = attribut' * pinv(PhiBPhi_new + lambda * eye(sampleTimeLength));
        
        Hy_new = H_new * label;
        Hy_new = Hy_new';
        Hy_new = reshape(Hy_new, m*nodeNum*usedNum, 1);
        
        HPhi_new = H_new*attribut;
        
        mu_x = zeros(nodeNum*usedNum*m,1);
        Sigma_x_new = zeros(usedNum*m,usedNum*m);
        Cov_x_new = [];
        
        currentLoc = 0;
        for i = 1 : usedNum
            currentLen = size(Sigma0{i},1);
            currentLoc = currentLoc + 1;
            seg = currentLoc : 1 : currentLoc + currentLen - 1;
            real_index = m*(i-1)+1;
            mu_x(seg) = Sigma0_new(real_index,real_index) * Hy_new(seg);% solution
            
            real_index_start = real_index;
            real_index_end = m*i;
            Cov_x_tmp = [];
            for index_inner = real_index_start : real_index_end
                %Sigma_x_new(i,i) = Sigma0_new(i,i) - Sigma0_new(i,i) * HPhi_new(i,i)*Sigma0_new(i,i); % solution Here need to fix
                Sigma_x_new(index_inner,index_inner) = Sigma0_new(index_inner,index_inner) - Sigma0_new(index_inner,index_inner) * HPhi_new(index_inner,index_inner)*Sigma0_new(index_inner,index_inner); % solution Here need to fix
                sigma_x_tmp = Sigma_x_new(index_inner,index_inner)*ones(1,currentLen/m);
                Cov_x_tmp = [Cov_x_tmp,sigma_x_tmp];
            end
            %Cov_x_new{i} = Sigma_x_new(i,i)*eye(currentLen,currentLen) + mu_x(seg) * mu_x(seg)';
            Cov_x_new{i} = diag(Cov_x_tmp) + mu_x(seg) * mu_x(seg)';
            currentLoc = seg(end);
        end
        
        
        % estimate gamma(i) and lambda
        gamma_old = gamma;
        lambdaComp = 0;
        for i =  1 : usedNum
            gamma(i) = trace(Cov_x_new{i})/size(Cov_x_new{i},1);
        end
        gamma_embd = [];
        for k = 1:usedNum
            for kk = 1 : m
                gamma_embd = [gamma_embd,gamma(k)];
            end
        end
        Sigma0_new = diag(gamma_embd);
        lambdaComp_new = trace(Sigma_x_new / sparse(Sigma0_new));
        mu_x_trans = reshape(mu_x,nodeNum,usedNum*m );
        mu_x_trans = mu_x_trans';
        tmp = attribut * mu_x_trans;
        lambda_new = norm(label - attribut * mu_x_trans,2)^2/TN + lambda_new * (NN - lambdaComp_new)/TN;
        lambda = lambda_new;
        
        
        
        % ================= Check stopping conditions, eyc. ==============
        if (size(mu_x) == size(mu_old))
            dmu = max(max(abs(mu_old - mu_x)));
            if (dmu < EPSILON)  break;  end;
        end;
        if (PRINT)
            disp([' iters: ',num2str(count),...
                ' num coeffs: ',num2str(usedNum), ...
                ' min gamma: ', num2str(min(gamma)),...
                ' gamma change: ',num2str(max(abs(gamma - gamma_old))),...
                ' mu change: ', num2str(dmu)]);
        end;
        if (count >= MAX_ITERS), if PRINT, fprintf('Reach max iterations. Stop\n\n'); end; break;  end;
        %RunTi{2} = [RunTi{2},toc];
    end;
end;

%% Expand hyperparameyers
gamma_used = sort(keep_list);
gamma_est = zeros(nodeNum,1);
gamma_est(keep_list,1) = gamma;


%% reconstruct the original signal
x = zeros(NN,1);
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
Result.gamma_est = gamma_est;
Result.count = count;
Result.lambda = lambda;
return;

