function [ attribute_comp,label_comp,transNet] = f_dynamicsGen_compac_discretet(sampleTimeLength,  networks,periodInterval, BER)
rest = mod(sampleTimeLength, periodInterval);
cycleNum = floor(sampleTimeLength/periodInterval);
[nodeNum,~] = size(networks);

attributs = [];
label = [];

tempAttr = zeros(periodInterval+1,nodeNum);


for i =1:(cycleNum)
    tempAttr(1,:) = rand(1,nodeNum);
    tempAttr(1,:) = f_trans2_01(tempAttr(1,:));
    for j = 1:periodInterval
        h = -1 *  tempAttr(j,:) * networks;
        sig_h=1./(1+exp(h));
        tempAttr(j+1,:) = f_trans2_01(sig_h);
    end
    attributs = [attributs; tempAttr(1:periodInterval,:)];
    templabel =  tempAttr(2:periodInterval+1,:);
    label = [label;templabel ];
end

if rest ~= 0
    tempAttr(1,:) = rand(1,nodeNum);
    tempAttr(1,:) = f_trans2_01(tempAttr(1,:));
    for j = 1:rest
        h = -1 *  tempAttr(j,:) * networks;
        sig_h=1./(1+exp(h));
        tempAttr(j+1,:) = f_trans2_01(sig_h);
    end
    
    attributs = [attributs; tempAttr(1:rest,:)];
    templabel =  tempAttr(2:rest+1,:);
    label = [label;templabel ];
    
end

for i = 1 : sampleTimeLength
    for j = 1 : nodeNum
        if rand<BER
            label(i,j) = abs(1-label(i,j));
        end
    end
end

%%%%%%%%%%% outputs %%%%%%%%%%%%%%
attribute_comp = attributs;
label_comp = label;
networks = networks';
transNet = reshape(networks, nodeNum*nodeNum,1);

