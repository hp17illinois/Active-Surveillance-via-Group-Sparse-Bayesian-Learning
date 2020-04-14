function SensorPosition = f_MPME_greed_forb(MeasurementMatrix,budget_k,forb_list)
%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
SensorPosition = zeros(budget_k,1);

%% ************* Determine the first to n-th sensor position ********************
    
%************* Determine the 1st sensor position ********************
RowNorm = sum(MeasurementMatrix.^2,2);
RowNorm(forb_list) = 0;
[~, SensorPosition(1)]=max(RowNorm);
%*********** determine 2nd to n-th sensor position ****************
Phi= MeasurementMatrix(SensorPosition(1),:);
OrthBasis = orth(Phi');
P = eye(n,n)-OrthBasis*OrthBasis';
MeasurementMatrix(SensorPosition(1),:)=zeros(1,n); 
for i = 2:budget_k
    ProjectionTemp = P*MeasurementMatrix';
    ProjectionTempNorm = sum(ProjectionTemp.^2);
    ProjectionTempNorm(forb_list)=0;
    [~, SensorPosition(i)]=max(ProjectionTempNorm);  
    
    Phi=[Phi;MeasurementMatrix(SensorPosition(i),:)];
    OrthBasis = orth(Phi');
    P = eye(n,n)-OrthBasis*OrthBasis';    
    MeasurementMatrix(SensorPosition(i),:)=zeros(1,n);
end
