function SensorPosition = MPME(MeasurementMatrix,budget_k)

%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
% if (n > N); 
%  error('Rows is not enough'); 
% end
% if (budget_k < n); 
%  error('More sensors are needed'); 
% end

%% ******************* Initilization *******************************
SensorPosition = zeros(budget_k,1);

%% ************* Determine the first to n-th sensor position ********************
    
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
[Maxvalue, SensorPosition(1)]=max(RowNorm);
%*********** determine 2nd to n-th sensor position ****************
Phi= MeasurementMatrix(SensorPosition(1),:);
OrthBasis = orth(Phi');
P = eye(n,n)-OrthBasis*OrthBasis';
MeasurementMatrix(SensorPosition(1),:)=zeros(1,n); 
for i = 2:budget_k
    ProjectionTemp = P*MeasurementMatrix';
    ProjectionTempNorm = sum(ProjectionTemp.^2);
    [Maxvalue, SensorPosition(i)]=max(ProjectionTempNorm);  
    
    Phi=[Phi;MeasurementMatrix(SensorPosition(i),:)];
    OrthBasis = orth(Phi');
    P = eye(n,n)-OrthBasis*OrthBasis';    
    MeasurementMatrix(SensorPosition(i),:)=zeros(1,n);
end
