clear all;
clc

sig = 10; %Variance of gaussian noise that is artificially added
k = 3;
patchsize = 5;

%Reading the clean image
im_original = read_and_disp('images/barco.png');

%Adding artificial noise
[im_noisy, gaussian_noise] = add_gaussian_noise(im_original, sig);

%Wavelet decompostion
%Daubechies Wavelet decomposition
[cA1,cH1,cV1,cD1] = dwt2(im_noisy,'db2');
subbands = cat(3,cA1,cH1,cV1,cD1);
num_of_subbands = size(subbands);



%iterating over subbands
for i= 2:num_of_subbands(3)
%     res_subbands(:,:,i) = denoising_utility(subbands(:,:,i),k,patchsize,gaussian_noise); 
    res_subbands(:,:,i) = denoi_BLS_GSM_band(subbands(:,:,i),[5 5],gaussian_noise,0,1,1,10);
end
cA1_res = cA1;

cH1_res = res_subbands(:,:,2);
cV1_res = res_subbands(:,:,3);
cD1_res = res_subbands(:,:,4);
res_image = idwt2(cA1_res,cH1_res,cV1_res,cD1_res,'db2');
figure, imshow(int8(cH1));
figure, imshow(int8(cH1_res));
figure, imshow(int8(cV1));
figure, imshow(int8(cV1_res));
figure, imshow(int8(cD1));
figure, imshow(int8(cD1_res));
figure, imshow(res_image);
