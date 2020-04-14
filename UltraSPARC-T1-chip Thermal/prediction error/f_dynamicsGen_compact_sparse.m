function [Phi,Y,attribute_cons,label_cons,transNet] = f_dynamicsGen_compact_sparse(sampleTimeLength, networks, periodInterval)
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
label_cons = label/scaler;
%%%%%%%%%%% outputs %%%%%%%%%%%%%%
attribute_cons = attributs;

I = sparse(diag(ones(nodeNum,1)));
Phi = kron(attributs, I);
Phi = sparse(Phi);

signal = label_cons';
Y = reshape(signal, nodeNum * sampleTimeLength, 1);

networks = networks';
transNet = reshape(networks, nodeNum*nodeNum,1);

