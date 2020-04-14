function SensorPosition = MNEP(MeasurementMatrix,budget_k)
%% ************************  Introduction ********************************
% MNEP:              Minimum Nonzero Eigenvalue Pursuit
%
% For detail about this algorithm, See:
%
%     C. Jiang, Y. C. Soh and H. Li, Sensor placement by maximal projection on
%     minimum eigenspace for linear inverse problems, 2015. (submitted to IEEE-
%     TSP).
%
% INPUT:
% MeasurementMatrix: a matrix from which we need to choose 'sensorNO' rows to
%                    construct a new matrix, and the chosen row index corresponds 
%                    to the sensor positions which are saved in 'SensorPosition'.
% sensorNO:          sensor number
%
%
% OUTPUT:
% SensorPosition:    a vector whose elements are the indices of the rows
%                    of 'MeasurementMatrix' chosen to place sensor
%
%
%
% AUTHOR Information:
% Jiang Chaoyang, EEE, NTU 
% Email: cjiang003@e.ntu.edu.sg, chaoyangjiang@hotmail.com
% Finished at 10-sep-zo14 
% Modified at 15-Nov-2015
%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
% if (n > N); 
%  error('Rows is not enough'); 
% end
% if (budget_k < n); 
%  error('More sensors are needed'); 
% end

SensorPosition = zeros(budget_k,1);

%% ************* Determine the first sensor position ********************
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
[Maxvalue, SensorPosition(1)]=max(RowNorm);

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
        
        [Maxvalue, MaxvalueIndex]=max(s(1,:));
        SensorPosition(i) = s(2,MaxvalueIndex);        
        Phi(i,:)= MeasurementMatrix(SensorPosition(i),:);          
    end   
end