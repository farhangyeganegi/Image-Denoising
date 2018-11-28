clc
clear
close all 
%%%%%%%%%%%%%%%%%%%%%%%%%%
I = imread('90.jpg');
l = size(I);
A = zeros(l);
for i=1:l(1)
    for j=1:l(2)
        A(i,j) = I(i,j);
    end
end
sigma = input('Please enter the value of the standard deviation of the gaussian noise = ');
n = sigma*randn(l);
An = A + n;
% applying SVD to the noisy image;
[s1, v1, d1] = svd(An);
T = diag(v1)';
beta = 1:0.5:60;
M = zeros(1,length(beta));
p = 1;
for beta=1:0.5:60
    for i=length(T):-1:1
    if T(i)>=beta*sigma
        e = i;
        break
    end
    end
    s11 = s1(1:l(1),1:e);
    v11 = v1(1:e,1:e);
    k1 = d1';
    d11 = k1(1:e,1:l(2));
    A_reconstruct = s11*v11*d11;
    M(p) = MSE(A_reconstruct,A);
    p = p + 1;
end
beta = 1:0.5:60;
figure;
plot(beta,M,'*');
xlabel('beta');
ylabel('MSE');
title('mean-square-error');

