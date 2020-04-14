function SensorPosition = f_MNEP_greed_modified(MeasurementMatrix,budget_k)
[N,n] = size(MeasurementMatrix);  
del_num= N - budget_k;
del_list = [1:N];
% if (n > N); 
%  error('Rows is not enough'); 
% end
% if (budget_k < n); 
%  error('More sensors are needed'); 
% end

SensorPosition = zeros(del_num,1);

%% ************* Determine the first sensor position ********************
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
[Maxvalue, SensorPosition(1)]=max(RowNorm);
del_index = del_list == SensorPosition(1);
del_list(del_index)=[];

%% ********* determine the senond to nth sensor positions
Phi = zeros(budget_k,n);
Phi(1,:)=MeasurementMatrix(SensorPosition(1),:);
if budget_k>1
    for i = 2:del_num 
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
        
        [Maxvalue, MaxvalueIndex]=max(s(1,:));
        SensorPosition(i) = s(2,MaxvalueIndex);
        del_index = del_list == SensorPosition(i);
        del_list(del_index)=[];
        Phi(i,:)= MeasurementMatrix(SensorPosition(i),:);
    end
    SensorPosition = del_list;
else
    SensorPosition = SensorPosition(1);
end
