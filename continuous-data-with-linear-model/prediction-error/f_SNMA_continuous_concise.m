function Result = f_SNMA_continuous_concise(attribut,label, nodeNum,sampleTimeLength,sen_k)
%global RunTi;
% Default Parameter Values for Any Cases
EPSILON       = 1e-10;       % solution accurancy tolerance
MAX_ITERS     = 500;        % maximum iterations
lambda = 1e-3;
lambda_new = 1e-3;
PRINT =  0;

%% Initialization
TN = sampleTimeLength * nodeNum;
NN = nodeNum * nodeNum;

attribut0 = attribut;
blkLen = nodeNum;
blkStartLoc = [1:blkLen:NN];
blkStartLoc0 = blkStartLoc;
blkLenList = blkLen * ones(1,nodeNum);

maxLen = nodeNum;

Sigma0_new = eye(nodeNum);
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
        %=========== Prune weights as their \gamma go to zero ==============
        index = find(gamma >1e-2);
        usedNum = length(index);
        keep_list = keep_list(index);
        blkStartLoc = blkStartLoc(index);
        blkLenList = blkLenList(index);
        % prune gamma and associated nodes in Sigma0 and data
        gamma = gamma(index);
        temp = Sigma0;
        Sigma0_new = diag(gamma);
        Sigma0 = [];
        for k = 1 : usedNum
            Sigma0{k} = temp{index(k)};
        end
        
        % construct new attributes
        temp = [];
        for k = 1 : usedNum
            temp = [temp, attribut0(:,keep_list(k))];
        end
        attribut = temp;
    end
    while (1)

        %=================== Compute new weights =================
        mu_old = mu_x;
        
        PhiBPhi_new = attribut*Sigma0_new*attribut';
        H_new = attribut' /(PhiBPhi_new + lambda * eye(sampleTimeLength));
        
        Hy_new = H_new * label;
        Hy_new = Hy_new';
        Hy_new = reshape(Hy_new, nodeNum*usedNum, 1);
        
        HPhi_new = H_new*attribut;
        
        mu_x = zeros(nodeNum*usedNum,1);
        Sigma_x_new = zeros(usedNum,usedNum);
        Cov_x_new = [];
        
        invB = [];
        currentLoc = 0;
        for i = 1 : usedNum
            currentLen = nodeNum;
            currentLen = size(Sigma0{i},1);
            currentLoc = currentLoc + 1;
            seg = currentLoc : 1 : currentLoc + currentLen - 1;
            mu_x(seg) = Sigma0_new(i,i) * Hy_new(seg);% solution
            
            Sigma_x_new(i,i) = Sigma0_new(i,i) - Sigma0_new(i,i) * HPhi_new(i,i)*Sigma0_new(i,i);
            
            Cov_x_new{i} = Sigma_x_new(i,i)*eye(currentLen,currentLen) + mu_x(seg) * mu_x(seg)';
            currentLoc = seg(end);
            invB{i} = eye(currentLen);
        end
        
        
        % estimate gamma(i) and lambda
        
        gamma_old = gamma;
        lambdaComp = 0;
        for i =  1 : usedNum
            gamma(i) = trace(invB{i} * Cov_x_new{i})/size(Cov_x_new{i},1);
        end
        
        Sigma0_new = diag(gamma);
        inv_Sigma0_new = diag(1./gamma);
        lambdaComp_new = trace(Sigma_x_new * sparse(inv_Sigma0_new));
        mu_x_trans = reshape(mu_x,nodeNum,usedNum );
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
%         toc;
        %RunTi{2} = [RunTi{2},toc];
    end;
end;

gamma_used = sort(keep_list);
gamma_est = zeros(nodeNum,1);
gamma_est(keep_list,1) = gamma;

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

