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
im_noisy = im_noisy2;

daub_order = 2;
[pyr,pind] = buildWUpyr(im_noisy,Nsc,daub_order);

%Extract the wavelets
[numpyr, col] = size(pind);
start=1;
subbands = [];
for i=2:numpyr
    subband = pyr(start:(pind(i,1)*pind(i,2)));
    subband = reshape(subband, pind(i,1), pind(i,2));
    start = start+ (pind(i,1)*pind(i,2));
end
daub_order = 2;
[pyr,pind] = buildWUpyr(im_noisy,Nsc,daub_order);
%Wavelet decompostion
%Daubechies Wavelet decomposition
% [cA1,cH1,cV1,cD1] = dwt2(im_noisy,'db2');
% subbands = cat(3,cA1,cH1,cV1,cD1);
% num_of_subbands = size(subbands);





%iterating over subbands
for i= 1:numpyr-1;
    [p_k_ym_subband(:,:,i), Cov_k_subband(:,:,:,i), p_z_k(:,:,i)] = expec_maxim(subbands(:,:,i),k,patchsize,gaussian_noise);  
end

