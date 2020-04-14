function [a_set_new,c_set_new] = f_GPMI_submodual(Sigma,u_set,a_set,c_set)

%global RunTi;

c_size = length(c_set);
au_set = horzcat(a_set,u_set);

min_en = 10000000000000000000000000000;
select_node = -1;
for i=1:c_size
  y = c_set(i);
  yau_set = horzcat(y,au_set);
  ya_set = horzcat(y,a_set);

  Sigma_yau = Sigma(yau_set,yau_set);
  Sigma_ya = Sigma(ya_set,ya_set);
  
%   test_ = Sigma(au_set,au_set);
%   Gaussian_entropy(test_)
  
  H_yau = f_GPMI_Gaussian_entropy(Sigma_yau);
  H_ya  = f_GPMI_Gaussian_entropy(Sigma_ya);
  
  H_condi = H_yau - H_ya;
  
  if H_condi < min_en
      min_en = H_condi;
      select_node = y;
  end
end
if select_node ~= -1
    c_set_new = c_set(c_set~=select_node);
    a_set_new = horzcat(select_node,a_set);
else
    c_set_new = c_set;
    a_set_new = a_set;
end

