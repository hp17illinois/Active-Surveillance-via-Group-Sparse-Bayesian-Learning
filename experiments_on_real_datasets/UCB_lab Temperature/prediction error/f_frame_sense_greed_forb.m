function SensorPosition = f_frame_sense_greed_forb(MeasurementMatrix,sensorNO,forb_list)
%% ********************* Preparation *****************************
[N,n] = size(MeasurementMatrix);  
%% ************* Initialization ***********************************
Nindex = 1:1:N;
Lindex = zeros(1,N-sensorNO); % The index of removed rows

%% ****************** Calculate the square Matrix ***************** 
Psi2 = (MeasurementMatrix*MeasurementMatrix').^2; % Psi2 = (Phi*Phi^T).^2
%% ******************** First two deleted sensor locations *************
[~,forb_sen_num] = size(forb_list);
for i = 1: forb_sen_num
    Lindex(i)=forb_list(i);
    Psi2(Lindex(i),:)=zeros(1,N);
    Psi2(:,Lindex(i))=zeros(N,1);
end

[C,IndexTemp] = max(Psi2); 
[Maxvalue, Lindex(forb_sen_num+1)] = max(C);
Lindex(forb_sen_num+2) = IndexTemp (Lindex(forb_sen_num+1));
Psi2(Lindex(forb_sen_num+1),:)=zeros(1,N);
Psi2(:,Lindex(forb_sen_num+1))=zeros(N,1);
if Lindex(forb_sen_num+1) ~= Lindex(forb_sen_num+2)
    Psi2(Lindex(forb_sen_num+2),:)=zeros(1,N);
    Psi2(:,Lindex(forb_sen_num+2))=zeros(N,1);
    sta_idx = forb_sen_num+3;
else
    sta_idx = forb_sen_num+2;
end
%% ** Determine Lindex (Remove the other rows one-by-one)*****************
for i = sta_idx:N-sensorNO
    OneRowInfluence = 2*sum(Psi2,2)-diag(Psi2);
    [Maxvalue, maxindex] = max(OneRowInfluence);
    Lindex(i)=maxindex;    

    Psi2(Lindex(i),:)=zeros(1,N);
    Psi2(:,Lindex(i))=zeros(N,1);
end
SensorPosition = setdiff(Nindex,Lindex);
