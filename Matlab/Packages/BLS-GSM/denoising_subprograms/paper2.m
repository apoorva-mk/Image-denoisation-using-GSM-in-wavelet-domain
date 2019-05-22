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



