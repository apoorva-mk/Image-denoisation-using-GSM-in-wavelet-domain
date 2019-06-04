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
%peaksnr = psnr(im0,im);


%dwt transform
%first vectorise
%im01 = im0(:);

%now apply wavelet transform
%[cA,cD] = dwt(im01,'sym4');

%cA and cD are the low and high pass filters
%X = idwt(cA,cD,'sym4');

%reshape
%x = reshape(X,[l,h]);
%x = uint8(x);
%figure,imshow(x);

%to view individual components cA and cD, need to upsample
%ca = upsample(cA,2);
% add 2-1 zeroes intermediately
% might have to delete some values
%then reshape and display
%ca = ca(1:l*h);
%y_new = reshape(ca,[l,h]);
%figure,imshow(uint8(y_new));



%Daubechies Wavelet decomposition
[cA1,cH1,cV1,cD1] = dwt2(im0,'db2');

%EM Implementation
patchsize = 5;

%padding with zeros to make it divisible by patch size
lpad = patchsize - mod(l,patchsize);
hpad = patchsize - mod(h,patchsize);
paddedim = padarray(im, [lpad,hpad], 'post'); 

%Vectorising the neighbourhoods
%Extracting patches and vectorising the patches and making a new matrix
patchmatrix = [];
for i = 1:patchsize:h-patchsize+1
    for j = 1:patchsize:l-patchsize+1
        h_end = i+patchsize-1;
        l_end = j+patchsize-1;
        patch = paddedim(i:h_end , j:l_end);
        patch = reshape(patch,[],1);
        patchmatrix = [ patchmatrix patch];
    end
end

%m stores the number of non-overlapping patches
[ patchlen, m ] = size(patchmatrix);

%Initialisation of proabability vectors 
k = 3;
numz = 13;
P_k_n = ones(1,k)/k;
p_z_k_n = ones(numz,k)/numz;

%Initializing z values
lzmin = -20.5;
lzmax = 3.5;
step = 2;
lzi = lzmin:step:lzmax;
numz = length(lzi);
zi = exp(lzi);
p_z = ones(1,numz)/numz;

%Next value vectors intialisation
P_k_n1 = zeros(1,k);
p_z_k_n1 = zeros(numz,k)/numz;


%Initialising the covariance matrices  !!!!!!!Change the intialisation!!!!!
Cov_k = [];
covcount = 0;
for i = 1:patchsize:h-patchsize+1
    for j = 1:patchsize:l-patchsize+1
        h_end = i+patchsize-1;
        l_end = j+patchsize-1;
        patch = paddedim(i:h_end , j:l_end);
        covtemp = cov(patch);
        Cov_k(:, :, covcount+1) = covtemp ;
        covcount = covcount + 1;
        if covcount == k
            break
        end
    end
    if covcount == k
        break
    end
end   


%EM Algorithm implementation
%Estimating the parameters

while( abs(sum(P_k_n)-sum(P_k_n1)) > 0.001)
    
%Calculating other probabilities that are required
%Calculating p(y(m)/k)

p_ym_k = zeroes(m, k);
p_ym = zeroes(1,m);
for i = 1:m
    for j = 1:k
        p_ym_k(i,j) = prob_y_given_k( patchmatrix(1:patchlen, i), Cov_k(:,:,j), p_z);
    end
    
    p_ym(i) = prob_y( p_ym_k( i, :), P_k_n);
end

%Equation 1:
%P_k_n1 = times (P_k_n , sum(p_y_k_n)) /  ( ( P_k_n * sum(p_y_k_n)) * m);
P_k_n1 =  times ( (P_k_n).', ( sum(p_ym_k)./ p_ym )) / m; 

%Equation 2:
p_k_z_n1 = p
end




















