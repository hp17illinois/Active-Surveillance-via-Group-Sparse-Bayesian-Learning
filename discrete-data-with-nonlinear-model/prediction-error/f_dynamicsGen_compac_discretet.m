function [ attribute_comp,label_comp,transNet] = f_dynamicsGen_compac_discretet(sampleTimeLength,  networks, BER)
[nodeNum,duration] = size(networks);

attributs = [];
label = [];

for i =1:(sampleTimeLength)
    tempAttr = rand(1,nodeNum);
    tempAttr = f_trans2_01(tempAttr);
    h = -1 * tempAttr * networks;
    sig_h=1./(1+exp(h));
    templabel = f_trans2_01(sig_h);
    attributs = [attributs; tempAttr];
    label = [label;templabel ];
end

for i = 1 : nodeNum
    for j = 1 : duration
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

