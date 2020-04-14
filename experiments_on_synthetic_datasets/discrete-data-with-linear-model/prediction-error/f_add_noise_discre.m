function [ y ] = f_add_noise_discre( ori_Y, BER)
y = ori_Y;
[n,m] = size(ori_Y);
if m == 1
    for i=1:n
        if rand()<BER
            y(i) = abs(1-ori_Y(i));
        end
    end
else
    for i = 1:n
        for j = 1:m
            if rand()<BER
                y(i,j) = abs(1-ori_Y(i,j));
            end
        end
    end
end

