function [ attributs,attributsEmbd, label , transNet] = f_dynamicsGen_compact_poly_fun_embd(sampleTimeLength, networks,m)
[nodeNum,~] = size(networks);
m=m+1;
nodeNum = nodeNum/m;
attributs = [];
attributsEmbd = [];
label = [];

% tempAttr = zeros(1,nodeNum);
% tempEmbd_1 = zeros(1,nodeNum);

for i =1:(sampleTimeLength)
    tempTotalAttr = [];
    tempAttr = rand(1,nodeNum);
    tempEmbd_1 = 0.5*(tempAttr.^2); % basic function is x and x^2
    for k =1:nodeNum
        tempTotalAttr = [tempTotalAttr, tempAttr(k) ];
        tempTotalAttr = [tempTotalAttr, tempEmbd_1(k) ];
    end
    templabel = tempTotalAttr * networks;
    attributs = [attributs; tempAttr];
    attributsEmbd = [attributsEmbd; tempTotalAttr];
    label = [label;templabel ];
end

scaler = max(max(attributs));
attributs = attributs/scaler;
attributsEmbd = attributsEmbd/scaler;
label = label/scaler;

networks = networks';
transNet = reshape(networks, m*nodeNum*nodeNum,1);

