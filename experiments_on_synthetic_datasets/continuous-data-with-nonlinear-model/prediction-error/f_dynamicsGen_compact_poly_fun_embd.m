function [ attributs,attributsEmbd, label , transNet] = f_dynamicsGen_compact_poly_fun_embd(sampleTimeLength, networks,m,periodInterval)

rest = mod(sampleTimeLength, periodInterval);
cycleNum = floor(sampleTimeLength/periodInterval);
[nodeNum,~] = size(networks);

m=m+1;
nodeNum = nodeNum/m;
attributs = [];
attributsEmbd = [];
label = [];

% tempAttr = zeros(1,nodeNum);
% tempEmbd_1 = zeros(1,nodeNum);
tempAttr = zeros(periodInterval+1,nodeNum);
tempEmbd_1 = zeros(periodInterval+1,m*nodeNum);
%attr_index = 1:m:m*nodeNum;
for i =1:(cycleNum)
    tempAttr(1,:) = rand(1,nodeNum);
    
    for j = 1:periodInterval
        tempEmbd = 2*(tempAttr(j,:).^2); % basic function is x and x^2
        tempTotalAttr = [];
        for k =1:nodeNum
            tempTotalAttr = [tempTotalAttr, tempAttr(j,k) ];
            tempTotalAttr = [tempTotalAttr, tempEmbd(k) ];
        end
        tempEmbd_1(j,:) = tempTotalAttr;
        tempAttr(j+1,:) = tempEmbd_1(j,:) * networks;
    end
    
    attributsEmbd = [attributsEmbd; tempEmbd_1(1:periodInterval,:)];
    attributs = [attributs; tempAttr(1:periodInterval,:)];
    templabel =  tempAttr(2:periodInterval+1,:);
    label = [label;templabel];
end

scaler = max(max(attributs));
attributs = attributs/scaler;
attributsEmbd = attributsEmbd/scaler;
label = label/scaler;

networks = networks';
transNet = reshape(networks, m*nodeNum*nodeNum,1);

