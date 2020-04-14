function [ y ] = f_add_noise_continu( ori_Y, SNR)

[n,m] = size(ori_Y);
if m == 1
    stdnoise = std(ori_Y)*10^(-SNR/20);
    noise = randn(n,1) * stdnoise;
    y = ori_Y + noise;
else
    ori_Y = ori_Y';
    ori_Y = reshape(ori_Y, n * m, 1);
    stdnoise = std(ori_Y)*10^(-SNR/20);
    noise = randn(n*m,1) * stdnoise;
    y = ori_Y + noise;
    y = reshape(y,m,n );
    y = y';
end

end

