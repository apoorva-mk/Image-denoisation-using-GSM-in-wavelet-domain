%Generating noisy image
function [im_noisy,noise] = add_gaussian_noise(im0, sig) 

seed = 0;
randn('state', seed);
noise = randn(size(im0));
noise = noise/sqrt(mean2(noise.^2));
im_noisy = im0 + sig*noise; 

figure(2)
rang = showIm(im_noisy,'auto');title('Noisy Image');
