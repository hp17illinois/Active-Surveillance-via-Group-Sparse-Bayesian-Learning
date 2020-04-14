function J = f_group_lasso_discrete_costFunction( inf_w,x,y,A )
[n,~]=size(y);
J = 0;
w=inf_w;
h=1./(1+exp(-x*inf_w));

J = (h-y)'*(h-y) + w'*A*w;

end


