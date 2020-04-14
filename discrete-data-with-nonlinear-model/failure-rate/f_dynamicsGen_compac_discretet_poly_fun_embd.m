function [ attributs,attributsEmbd, label , transNet] = f_dynamicsGen_compac_discretet_poly_fun_embd(sampleTimeLength, networks,m,BER)
[nodeNum,~] = size(networks);
m=m+1;
nodeNum = nodeNum/m;

attributs = [];
attributsEmbd = [];
label = [];

for i =1:(sampleTimeLength)
    tempTotalAttr = [];
    tempAttr = rand(1,nodeNum);
    tempAttr = f_trans2_01(tempAttr);
    tempEmbd_1 = 0.5*(tempAttr.^2); % basic function is x and x^2
    for k =1:nodeNum
        tempTotalAttr = [tempTotalAttr, tempAttr(k) ];
        tempTotalAttr = [tempTotalAttr, tempEmbd_1(k) ];
    end
    h = -1 * tempTotalAttr * networks;
    sig_h=1./(1+exp(h));
    templabel = f_trans2_01(sig_h);
    
    attributs = [attributs; tempAttr];
    attributsEmbd = [attributsEmbd; tempTotalAttr];
    label = [label;templabel ];
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

