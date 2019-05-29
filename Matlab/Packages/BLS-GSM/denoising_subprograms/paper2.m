clear
im0 = imread('images/barco.png');
im0 = double(im0);
figure(1)
rang0 = showIm(im0,'auto');title('Original image');

%Generating noisy image
sig = 10;
seed = 0;
randn('state', seed);
noise = randn(size(im0));
noise = noise/sqrt(mean2(noise.^2));
im = im0 + sig*noise; 

figure(2)
rang = showIm(im,'auto');title('Noisy Image');

[l,h]= size(im0);

%calculating psnr
peaksnr = psnr(im0,im);


%dwt transform
%first vectorise
im01 = im0(:);

%now apply wavelet transform
[cA,cD] = dwt(im01,'sym4');

%cA and cD are the low and high pass filters
X = idwt(cA,cD,'sym4');

%reshape
x = reshape(X,[l,h]);
x = uint8(x);
figure,imshow(x);

%to view individual components cA and cD, need to upsample
ca = upsample(cA,2);
% add 2-1 zeroes intermediately
% might have to delete some values
%then reshape and display
ca = ca(1:l*h);
y_new = reshape(ca,[l,h]);
figure,imshow(uint8(y_new));



%Daubechies Wavelet decomposition
[cA1,cH1,cV1,cD1] = dwt2(im0,'db2');

%EM Implementation
%Initialisation 
k = 10;
numz = 13;
Pk = ones(1,k)/k;
pk = ones(k,numz)/(k*numz);
Cov









