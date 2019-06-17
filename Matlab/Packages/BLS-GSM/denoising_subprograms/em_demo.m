clear all;
clc

sig = 10; %Variance of gaussian noise that is artificially added
k = 10;
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
    [p_k_ym_subband(:,:,i), Cov_k_subband(:,:,:,i), p_z_k(:,:,i)] = expec_maxim(subbands(:,:,i),k,patchsize,gaussian_noise);  
end


