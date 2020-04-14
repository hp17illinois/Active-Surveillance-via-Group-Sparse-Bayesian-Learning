function SensorPosition = f_MNEP_greed_forb(MeasurementMatrix,budget_k,feasible_list)
%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
SensorPosition = zeros(budget_k,1);

%% ************* Determine the first sensor position ********************
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
RowNorm_feas = RowNorm(feasible_list);
[~, max_index]=max(RowNorm_feas);
SensorPosition(1) = max_index;
del_index = feasible_list==SensorPosition(1);
feasible_list(del_index)=[];
%% ********* determine the senond to nth sensor positions
Phi = zeros(budget_k,n);
Phi(1,:)=MeasurementMatrix(SensorPosition(1),:);
if budget_k>1
    for i = 2:budget_k 
        s = [];
        for j = 1:N
            flag =1;
            for k = 1:i-1
                if j==SensorPosition(k) %****
                    flag = 0;
                    break;
                end
            end
            
            if flag
                Phi(i,:)= MeasurementMatrix(j,:);
                if i<n
                    eigenvalue = abs(eig(Phi(1:i,:)*Phi(1:i,:)'));
                else
                    eigenvalue = abs(eig(Phi(1:i,:)'*Phi(1:i,:)));
                end
                s = [s [min(eigenvalue);j]];               
            end
        end
        [~,row_num] = size(s);
        s_feas = [];
        for ind = 1:row_num
            if ismember(s(2,ind),feasible_list)
                s_feas=  [s_feas,s(:,ind)];
            end
        end
        [~, MaxvalueIndex]=max(s_feas(1,:));
        SensorPosition(i) = s(2,MaxvalueIndex);
        del_index = feasible_list==SensorPosition(i);
        feasible_list(del_index)=[];
        Phi(i,:)= MeasurementMatrix(SensorPosition(i),:);
    end   
end