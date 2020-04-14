function [ attributs,attributsEmbd, label , transNet] = f_dynamicsGen_compac_discretet_poly_fun_embd(sampleTimeLength, networks,m,periodInterval,BER)
rest = mod(sampleTimeLength, periodInterval);
cycleNum = floor(sampleTimeLength/periodInterval);


[nodeNum,~] = size(networks);
m=m+1;
nodeNum = nodeNum/m;

attributs = [];
attributsEmbd = [];
label = [];
tempAttr = zeros(periodInterval+1,nodeNum);
tempEmbd = zeros(periodInterval+1,nodeNum);
tempTotalAttr =  zeros(periodInterval+1,m*nodeNum);
for i =1:(cycleNum)
    tempAttr(1,:) = rand(1,nodeNum)-0.4; % small number of cases in the begining
    tempAttr(1,:) = f_trans2_01(tempAttr(1,:));
    tempEmbd(1,:) = zeros(1,nodeNum); % basic function is x_t - x_(t-1)
    
    for j = 1:periodInterval
        tempTotalAttr_ = [];
        for k =1:nodeNum
            tempTotalAttr_ = [tempTotalAttr_, tempAttr(j,k) ];
            tempTotalAttr_ = [tempTotalAttr_, tempEmbd(j,k) ];
        end
        tempTotalAttr(j,:) = tempTotalAttr_;
        
        h = -1 *  tempTotalAttr(j,:) * networks;
        sig_h=1./(1+exp(h));
        tempAttr(j+1,:) = f_trans2_01(sig_h);
        tempEmbd(j+1,:) = tempAttr(j+1,:) - tempAttr(j,:); % basic function is x_t - x_(t-1)
    end
    
    attributs = [attributs; tempAttr(1:periodInterval,:)];
    attributsEmbd = [attributsEmbd;tempTotalAttr(1:periodInterval,:) ];
    label     = [label;     tempAttr(2:periodInterval+1,:)];
end

if rest ~= 0
    tempAttr(1,:) = rand(1,nodeNum)-0.4; % small number of cases in the begining
    tempAttr(1,:) = f_trans2_01(tempAttr(1,:));
    tempEmbd(1,:) = zeros(1,nodeNum); % basic function is x_t - x_(t-1)
    
    for j = 1:rest
        tempTotalAttr_ = [];
        for k =1:nodeNum
            tempTotalAttr_ = [tempTotalAttr_, tempAttr(j,k) ];
            tempTotalAttr_ = [tempTotalAttr_, tempEmbd(j,k) ];
        end
        tempTotalAttr(j,:) = tempTotalAttr_;
        
        h = -1 *  tempTotalAttr(j,:) * networks;
        sig_h=1./(1+exp(h));
        tempAttr(j+1,:) = f_trans2_01(sig_h);
        tempEmbd(j+1,:) = tempAttr(j+1,:) - tempAttr(j,:); % basic function is x_t - x_(t-1)
    end
    
    attributs = [attributs; tempAttr(1:rest,:)];
    attributsEmbd = [attributsEmbd,tempTotalAttr(1:rest,:) ];
    label     = [label;     tempAttr(2:rest+1,:)];
end


for i = 1 : sampleTimeLength
    for j = 1 :  nodeNum
        if rand<BER
            label(i,j) = abs(1-label(i,j));
        end
    end
end

networks = networks';
transNet = reshape(networks, m*nodeNum*nodeNum,1);

