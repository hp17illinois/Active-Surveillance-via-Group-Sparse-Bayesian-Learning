function [ attributs, label , transNet] = f_dynamicsGen_compact(sampleTimeLength, networks, periodInterval )
rest = mod(sampleTimeLength, periodInterval);
cycleNum = floor(sampleTimeLength/periodInterval);
[nodeNum,~] = size(networks);

attributs = [];
label = [];

tempAttr = zeros(periodInterval+1,nodeNum);

for i =1:(cycleNum)
    tempAttr(1,:) = rand(1,nodeNum);
    for j = 1:periodInterval
        tempAttr(j+1,:) = tempAttr(j,:) * networks;
    end
    templabel =  tempAttr(2:periodInterval+1,:);
    attributs = [attributs; tempAttr(1:periodInterval,:)];
    label = [label;templabel ];
end

if rest ~= 0
    tempAttr(1,:) = rand(1,nodeNum);
    for j = 1:rest
        tempAttr(j+1,:) = tempAttr(j,:) * networks;
    end
    templabel =  tempAttr(2:rest+1,:);
    attributs = [attributs; tempAttr(1:rest,:)];
    label = [label;templabel ];
end
scaler = max(max(attributs));
attributs = attributs/scaler;
label = label/scaler;

networks = networks';
transNet = reshape(networks, nodeNum*nodeNum,1);

