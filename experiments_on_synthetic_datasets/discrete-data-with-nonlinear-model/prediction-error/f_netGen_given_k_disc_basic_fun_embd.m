function [ networks,hub_list ] = f_netGen_given_k_disc_basic_fun_embd(nodeNum, k_gt, sparsity,m)
m = m+1;
hub_list = [];
networks = zeros(nodeNum*m,nodeNum);
rowArr = randperm(nodeNum);
richRowNum = k_gt;

for i=1:richRowNum
    hub_list = [hub_list,rowArr(i)];
    row_num = (rowArr(i)-1)*2 +1;
    for j = 1:nodeNum
        temp = rand();
        sign_tmp = rand();
        if temp > sparsity
            if sign_tmp > 0.5
                networks(row_num,j) = rand * 100;
            else
                networks(row_num,j) = -1* rand * 100;
            end
        else
            networks(row_num,j) = 0;
        end
        temp = rand();
        sign_tmp = rand();
        if temp > sparsity
            if sign_tmp > 0.5
                networks(row_num+1,j) = rand * 100;
            else
                networks(row_num+1,j) =-1* rand * 100;
            end
        else
            networks(row_num+1,j) = 0;
        end
    end
end

for i=richRowNum+1:nodeNum
    row_num = (rowArr(i)-1)*2 +1;
    for j = 1:nodeNum
        temp = rand();
        sign_tmp = rand();
        if temp > 1 - 0.15 * (1-sparsity)
            if sign_tmp >0.5
                networks(row_num,j) = rand * 10;
            else
                networks(row_num,j) =-1* rand * 10;
            end
        else
            networks(row_num,j) = 0;
        end
        temp = rand();
        sign_tmp = rand();
        if temp > 1 - 0.15 * (1-sparsity)
            if sign_tmp >0.5
                networks(row_num+1,j) = rand * 10;
            else
                networks(row_num+1,j) = -1*rand * 10;
            end
        else
            networks(row_num+1,j) = 0;
        end
    end
end
ex = roundn(networks,-4);
networks = ex/max(max(ex));
hub_list = sort(hub_list);
% ex = textread('network.txt');
end