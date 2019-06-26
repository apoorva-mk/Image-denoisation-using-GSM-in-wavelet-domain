clear all;
clc

sig = 10; %Variance of gaussian noise that is artificially added
k = 10;
patchsize = 5;

%Reading the clean image
im_original = read_and_disp('images/barbara.png');



%Adding artificial noise
[im_noisy, gaussian_noise] = add_gaussian_noise(im_original, sig);

[Ny,Nx] = size(im_noisy);
Nsc = ceil(log2(min(Ny,Nx)) - 4);
hpad = (2^(Nsc+1)) - mod(Nx,(2^(Nsc+1)));
lpad = (2^(Nsc+1)) - mod(Ny, (2^(Nsc+1)));

im_noisy2 = padarray(im_noisy, [lpad,hpad], 'post'); 
gaussian_noise = padarray(gaussian_noise, [lpad,hpad], 'post');
im_noisy = im_noisy2;

daub_order = 2;
[pyr,pind] = buildWUpyr(im_noisy,Nsc,daub_order);

%Extract the wavelets
[numpyr, col] = size(pind);
start=1;
subbands = [];
% for i=2:numpyr
%     subband = pyr(start:(pind(i,1)*pind(i,2)));
%     subband = reshape(subband, pind(i,1), pind(i,2));
%     start = start+ (pind(i,1)*pind(i,2));
% end
% 
% %Wavelet decompostion
% %Daubechies Wavelet decomposition
% % [cA1,cH1,cV1,cD1] = dwt2(im_noisy,'db2');
% % subbands = cat(3,cA1,cH1,cV1,cD1);
% % num_of_subbands = size(subbands);
% 
% %iterating over subbands
% for i= 1:numpyr-1;
%     [p_k_ym_subband(:,:,i), Cov_k_subband(:,:,:,i), p_z_k(:,:,i)] = expec_maxim(subbands(:,:,i),k,patchsize,gaussian_noise);  
% end
i =2 ;
subband_1 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_1 = reshape(subband_1, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_2 = pyr(start: start+(pind(i,1)*pind(i,2))-1);
subband_2 = reshape(subband_2, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_3 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_3 = reshape(subband_3, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_4 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_4 = reshape(subband_4, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_5 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_5 = reshape(subband_5, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_6 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_6 = reshape(subband_6, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_7 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_7 = reshape(subband_7, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_8 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_8 = reshape(subband_8, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_9 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_9 = reshape(subband_9, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_10 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_10 = reshape(subband_10, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_11 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_11 = reshape(subband_11, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_12 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_12 = reshape(subband_12, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_13 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_13 = reshape(subband_13, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_14 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_14 = reshape(subband_14, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

subband_15 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_15 = reshape(subband_15, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;


subband_16 = pyr(start:start+(pind(i,1)*pind(i,2))-1);
subband_16 = reshape(subband_16, pind(i,1), pind(i,2));

start = start+(pind(i,1)*pind(i,2));
i = i+1;

%p_k_ym_subband_1 = denoising_utility(subband_1,k,patchsize,gaussian_noise); 
[p_k_ym_subband_1, Cov_k_subband_1, p_z_k_1] = expec_maxim(subband_1,k,patchsize,gaussian_noise); 
[p_k_ym_subband_2, Cov_k_subband_2, p_z_k_2] = expec_maxim(subband_2,k,patchsize,gaussian_noise); 
[p_k_ym_subband_3, Cov_k_subband_3, p_z_k_3] = expec_maxim(subband_3,k,patchsize,gaussian_noise); 
[p_k_ym_subband_4, Cov_k_subband_4, p_z_k_4] = expec_maxim(subband_4,k,patchsize,gaussian_noise);
[p_k_ym_subband_5, Cov_k_subband_5, p_z_k_5] = expec_maxim(subband_5,k,patchsize,gaussian_noise); 
[p_k_ym_subband_6, Cov_k_subband_6, p_z_k_6] = expec_maxim(subband_6,k,patchsize,gaussian_noise); 
[p_k_ym_subband_7, Cov_k_subband_7, p_z_k_7] = expec_maxim(subband_7,k,patchsize,gaussian_noise); 
[p_k_ym_subband_8, Cov_k_subband_8, p_z_k_8] = expec_maxim(subband_8,k,patchsize,gaussian_noise);
[p_k_ym_subband_9, Cov_k_subband_9, p_z_k_9] = expec_maxim(subband_9,k,patchsize,gaussian_noise); 
% % [p_k_ym_subband_10, Cov_k_subband_10, p_z_k_10] = expec_maxim(subband_10,k,patchsize,gaussian_noise); 
[p_k_ym_subband_11, Cov_k_subband_11, p_z_k_11] = expec_maxim(subband_11,k,patchsize,gaussian_noise); 
[p_k_ym_subband_12, Cov_k_subband_12, p_z_k_12] = expec_maxim(subband_12,k,patchsize,gaussian_noise);
% %[p_k_ym_subband_13, Cov_k_subband_13, p_z_k_13] = expec_maxim(subband_13,k,patchsize,gaussian_noise); 
% %[p_k_ym_subband_14, Cov_k_subband_14, p_z_k_14] = expec_maxim(subband_14,k,patchsize,gaussian_noise); 
% %[p_k_ym_subband_15, Cov_k_subband_15, p_z_k_15] = expec_maxim(subband_15,k,patchsize,gaussian_noise); 
% %[p_k_ym_subband_16, Cov_k_subband_16, p_z_k_16] = expec_maxim(subband_16,k,patchsize,gaussian_noise);
