function SensorPosition = f_MPME_greed_modified(MeasurementMatrix,budget_k)

%% ********************* Preparation *****************************
%MeasurementMatrix =MeasurementMatrix(:,1:10);
[N,n] = size(MeasurementMatrix); 
del_num= N - budget_k;
del_list = [1:N];
% if (n > N); 
%  error('Rows is not enough'); 
% end
% if (budget_k < n); 
%  error('More sensors are needed'); 
% end

%% ******************* Initilization *******************************
SensorPosition = zeros(del_num,1);

%% ************* Determine the first to n-th sensor position ********************
    
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
[Maxvalue, SensorPosition(1)]=max(RowNorm);
del_index = del_list == SensorPosition(1);
del_list(del_index)=[];
%*********** determine 2nd to n-th sensor position ****************
Phi= MeasurementMatrix(SensorPosition(1),:);
OrthBasis = orth(Phi');
P = eye(n,n)-OrthBasis*OrthBasis';
MeasurementMatrix(SensorPosition(1),:)=zeros(1,n); 
for i = 2:del_num
    ProjectionTemp = P*MeasurementMatrix';
    ProjectionTempNorm = sum(ProjectionTemp.^2);
    [Maxvalue, SensorPosition(i)]=max(ProjectionTempNorm);
    del_index = del_list == SensorPosition(i);
    del_list(del_index)=[];
    
    Phi=[Phi;MeasurementMatrix(SensorPosition(i),:)];
    OrthBasis = orth(Phi');
    P = eye(n,n)-OrthBasis*OrthBasis';    
    MeasurementMatrix(SensorPosition(i),:)=zeros(1,n);
end
SensorPosition = del_list';
