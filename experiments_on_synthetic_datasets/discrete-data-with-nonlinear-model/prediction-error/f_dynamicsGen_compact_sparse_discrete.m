function [Phi,Y,attribute_comp,label_comp,transNet] = f_dynamicsGen_compact_sparse_discrete(sampleTimeLength, networks, BER)
[nodeNum,duration] = size(networks);

attributs = [];
label = [];
tempAttr = zeros(1,nodeNum);

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

I = sparse(diag(ones(nodeNum,1)));
Phi = kron(attributs, I);
Phi = sparse(Phi);

label_comp = label;
signal = label';
Y = reshape(signal, nodeNum * sampleTimeLength, 1);

networks = networks';
transNet = reshape(networks, nodeNum*nodeNum,1);

