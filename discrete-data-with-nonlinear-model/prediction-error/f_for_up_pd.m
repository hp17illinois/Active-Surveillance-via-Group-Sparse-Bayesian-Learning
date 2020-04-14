function updat_w = f_for_up_pd(inf_w,x,y,A)

mu = 1./(1+exp(-x*inf_w));
pd_w = x' * (y-mu);
pd_w = pd_w - A*inf_w;

% V = mu.*(1-mu);
% V = diag(V);
% H = -x'*V*x-A;
% updat_w = H\pd_w;

updat_w = -pd_w;