function [ networks,hub_list ] = f_netGen_given_k(nodeNum, k_gt, sparsity)

hub_list = [];
networks = zeros(nodeNum);
rowArr = randperm(nodeNum);
richRowNum = k_gt;

for i=1:richRowNum
   hub_list = [hub_list,rowArr(i)];
   for j = 1:nodeNum
       temp = rand();
       if temp > sparsity
          networks(rowArr(i),j) = rand * 100;
       else
           networks(rowArr(i),j) = 0;
       end 
   end
end

for i=richRowNum+1:nodeNum
   for j = 1:nodeNum
       temp = rand();
       if temp > 1 - 0.15 * (1-sparsity)
          networks(rowArr(i),j) = rand * 10;
       else
           networks(rowArr(i),j) = 0;
       end 
   end
end
ex = roundn(networks,-4);
networks = ex/max(max(ex));
hub_list = sort(hub_list);
% ex = textread('network.txt');
end