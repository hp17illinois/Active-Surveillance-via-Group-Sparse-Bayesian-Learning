function [ out_vec ] = f_trans2_01( input_vec )
    input_vec(input_vec>0.5) = 1;
    input_vec(input_vec<=0.5) = 0;
    out_vec = input_vec;
end

