function  [attributs,attributsEmbd , label, duration, nodeNum] = loadRealData_lab_embd()

load_matrix = xlsread('Data_lab.csv','A1:AY203');
% load_matrix = xlsread('Data_lab_half_hour.csv','A1:AY408');
attributs = load_matrix(1:end-1,:);
label = load_matrix(2:end,:);
[duration, nodeNum] = size(attributs);

attributsEmbd = [];
for j = 1:duration
    tempEmbd = (attributs(j,:).^2); % basic function is x and x^2
    tempTotalAttr = [];
    for k =1:nodeNum
        tempTotalAttr = [tempTotalAttr, attributs(j,k) ];
        tempTotalAttr = [tempTotalAttr, tempEmbd(k) ];
    end
    attributsEmbd = [attributsEmbd;tempTotalAttr];
end
