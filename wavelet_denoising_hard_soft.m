% Image denoising using universal threshold in wavelet domain (hard thresholding);
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
sigma = 10; % noise level
B_1 = B_0 + sigma*randn(size(B_0));
wname = 'db10';
[CA,CH,CV,CD] = dwt2(B_1,wname);

p_1 = size(CA);
p_2 = size(CH);
p_3 = size(CV);
p_4 = size(CD);
CD_new = reshape(CD,1,p_4(1)*p_4(2));
sigma_HH = median(abs(CD_new))/0.6745;

T_H = sigma_HH*sqrt(2*log10(p_2(1)*p_2(2)));
T_V = sigma_HH*sqrt(2*log10(p_3(1)*p_3(2)));
T_D = sigma_HH*sqrt(2*log10(p_4(1)*p_4(2)));

for i=1:p_4(1)
    for j=1:p_4(2)
        if abs(CD(i,j))>=T_D
            CD(i,j) = CD(i,j);
        else
            CD(i,j) = 0;
        end
    end
end

for i=1:p_3(1)
    for j=1:p_3(2)
        if abs(CV(i,j))>=T_V
            CV(i,j) = CV(i,j);
        else
            CV(i,j) = 0;
        end
    end
end

for i=1:p_2(1)
    for j=1:p_2(2)
        if abs(CH(i,j))>=T_H
            CH(i,j) = CH(i,j);
        else
            CH(i,j) = 0;
        end
    end
end

B_2 = idwt2(CA,CH,CV,CD,'db10');
subplot(1,3,1);
imagesc(B_0);
colormap gray;
subplot(1,3,2);
imagesc(B_1);
colormap gray;
subplot(1,3,3);
imagesc(B_2);
colormap gray;
MSE_1=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_1 = MSE_1+((B_0(i,j)-B_1(i,j))^2)/(l(1)*l(2));
    end
end
MSE_2=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_2 = MSE_2+((B_0(i,j)-B_2(i,j))^2)/(l(1)*l(2));
    end
end
M = [MSE_1 MSE_2];
snr_1 = calSNR(B_0,B_1);
snr_2 = calSNR(B_0,B_2);
snr = [snr_1 snr_2];

%%
% Image denoising using universal threshold in wavelet domain (soft thresholding);
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
sigma = 10; % noise level
B_1 = B_0 + sigma*randn(size(B_0));
wname = 'db10';
[CA,CH,CV,CD] = dwt2(B_1,wname);

p_1 = size(CA);
p_2 = size(CH);
p_3 = size(CV);
p_4 = size(CD);
CD_new = reshape(CD,1,p_4(1)*p_4(2));
sigma_HH = median(abs(CD_new))/0.6745;

T_H = sigma_HH*sqrt(2*log10(p_2(1)*p_2(2)));
T_V = sigma_HH*sqrt(2*log10(p_3(1)*p_3(2)));
T_D = sigma_HH*sqrt(2*log10(p_4(1)*p_4(2)));

for i=1:p_4(1)
    for j=1:p_4(2)
        if CD(i,j)>=T_D
            CD(i,j) = CD(i,j)-T_D;
        elseif CD(i,j)<=-T_D
            CD(i,j) = CD(i,j)+T_D;
        else
            CD(i,j) = 0;
        end
    end
end

for i=1:p_3(1)
    for j=1:p_3(2)
        if CV(i,j)>=T_V
            CV(i,j) = CV(i,j)-T_V;
        elseif CV(i,j)<=-T_V
            CV(i,j) = CV(i,j)+T_V;
        else
            CV(i,j) = 0;
        end
    end
end

for i=1:p_2(1)
    for j=1:p_2(2)
        if CH(i,j)>=T_H
            CH(i,j) = CH(i,j)-T_H;
        elseif CH(i,j)<=-T_H
            CH(i,j) = CH(i,j)+T_H;
        else 
            CH(i,j) = 0;
        end
    end
end

B_2 = idwt2(CA,CH,CV,CD,'db10');
subplot(1,3,1);
imagesc(B_0);
colormap gray;
subplot(1,3,2);
imagesc(B_1);
colormap gray;
subplot(1,3,3);
imagesc(B_2);
colormap gray;
MSE_1=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_1 = MSE_1+((B_0(i,j)-B_1(i,j))^2)/(l(1)*l(2));
    end
end
MSE_2=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_2 = MSE_2+((B_0(i,j)-B_2(i,j))^2)/(l(1)*l(2));
    end
end
M = [MSE_1 MSE_2];
snr_1 = calSNR(B_0,B_1);
snr_2 = calSNR(B_0,B_2);
snr = [snr_1 snr_2];

%%
% Image despeckling using universal threshold in wavelet domain (hard thresholding);
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
B_1 = speck(B_0,0.5,3,3);
B_1 = log(B_1);
wname = 'db10';
[CA,CH,CV,CD] = dwt2(B_1,wname);

p_1 = size(CA);
p_2 = size(CH);
p_3 = size(CV);
p_4 = size(CD);
CD_new = reshape(CD,1,p_4(1)*p_4(2));
sigma_HH = median(abs(CD_new))/0.6745;

T_H = sigma_HH*sqrt(2*log10(p_2(1)*p_2(2)));
T_V = sigma_HH*sqrt(2*log10(p_3(1)*p_3(2)));
T_D = sigma_HH*sqrt(2*log10(p_4(1)*p_4(2)));

for i=1:p_4(1)
    for j=1:p_4(2)
        if abs(CD(i,j))>=T_D
            CD(i,j) = CD(i,j);
        else
            CD(i,j) = 0;
        end
    end
end

for i=1:p_3(1)
    for j=1:p_3(2)
        if abs(CV(i,j))>=T_V
            CV(i,j) = CV(i,j);
        else
            CV(i,j) = 0;
        end
    end
end

for i=1:p_2(1)
    for j=1:p_2(2)
        if abs(CH(i,j))>=T_H
            CH(i,j) = CH(i,j);
        else
            CH(i,j) = 0;
        end
    end
end

B_2 = idwt2(CA,CH,CV,CD,'db10');
B_2 = exp(B_2);
B_1 = exp(B_1);
subplot(1,3,1);
imagesc(B_0);
colormap gray;
subplot(1,3,2);
imagesc(B_1);
colormap gray;
subplot(1,3,3);
imagesc(B_2);
colormap gray;
% comparison with the noisy image;
MSE_1=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_1 = MSE_1+((B_0(i,j)-B_1(i,j))^2)/(l(1)*l(2));
    end
end
% comparison with the denoised image;
MSE_2=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_2 = MSE_2+((B_0(i,j)-B_2(i,j))^2)/(l(1)*l(2));
    end
end
M = [MSE_1 MSE_2];
snr_1 = calSNR(B_0,B_1);
snr_2 = calSNR(B_0,B_2);
snr = [snr_1 snr_2];

%%
% Image despeckling using universal threshold in wavelet domain (soft thresholding);
I = imread('90.jpg');
l = size(I);
B_0 = zeros(l(1),l(2));
for i=1:l(1)
    for j=1:l(2)
        B_0(i,j)=I(i,j);
    end
end
B_1 = speck(B_0,0.5,3,3);
B_1 = log(B_1);
wname = 'db10';
[CA,CH,CV,CD] = dwt2(B_1,wname);

p_1 = size(CA);
p_2 = size(CH);
p_3 = size(CV);
p_4 = size(CD);
CD_new = reshape(CD,1,p_4(1)*p_4(2));
sigma_HH = median(abs(CD_new))/0.6745;

T_H = sigma_HH*sqrt(2*log10(p_2(1)*p_2(2)));
T_V = sigma_HH*sqrt(2*log10(p_3(1)*p_3(2)));
T_D = sigma_HH*sqrt(2*log10(p_4(1)*p_4(2)));

for i=1:p_4(1)
    for j=1:p_4(2)
        if CD(i,j)>=T_D
            CD(i,j) = CD(i,j)-T_D;
        elseif CD(i,j)<=-T_D
            CD(i,j) = CD(i,j)+T_D;
        else
            CD(i,j) = 0;
        end
    end
end

for i=1:p_3(1)
    for j=1:p_3(2)
        if CV(i,j)>=T_V
            CV(i,j) = CV(i,j)-T_V;
        elseif CV(i,j)<=-T_V
            CV(i,j) = CV(i,j)+T_V;
        else
            CV(i,j) = 0;
        end
    end
end

for i=1:p_2(1)
    for j=1:p_2(2)
        if CH(i,j)>=T_H
            CH(i,j) = CH(i,j)-T_H;
        elseif CH(i,j)<=-T_H
            CH(i,j) = CH(i,j)+T_H;
        else 
            CH(i,j) = 0;
        end
    end
end

B_2 = idwt2(CA,CH,CV,CD,'db10');
B_2 = exp(B_2);
B_1 = exp(B_1);
subplot(1,3,1);
imagesc(B_0);
colormap gray;
subplot(1,3,2);
imagesc(B_1);
colormap gray;
subplot(1,3,3);
imagesc(B_2);
colormap gray;
% comparison with the noisy image;
MSE_1=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_1 = MSE_1+((B_0(i,j)-B_1(i,j))^2)/(l(1)*l(2));
    end
end
% comparison with the denoised image;
MSE_2=0;
for i=1:l(1)
    for j=1:l(2)
        MSE_2 = MSE_2+((B_0(i,j)-B_2(i,j))^2)/(l(1)*l(2));
    end
end
M = [MSE_1 MSE_2];
snr_1 = calSNR(B_0,B_1);
snr_2 = calSNR(B_0,B_2);
snr = [snr_1 snr_2];







