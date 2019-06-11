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
for i = 1:patchsize:(h-patchsize+1)
    for j = 1:patchsize:(l-patchsize+1)
        
        h_end = i+patchsize-1;
        l_end = j+patchsize-1;
        patch = paddedim(j:l_end , i:h_end);
        patch = reshape(patch,[],1);
        patchmatrix = [ patchmatrix patch];
    end
end

%m stores the number of non-overlapping patches
[ patchlen, m ] = size(patchmatrix);

%Finding C_w
[nv,nh,nb] = size(im0);
block = [ 5,5];
nblv = nv-block(1)+1;	% Discard the outer coefficients 
nblh = nh-block(2)+1;   % for the reference (centrral) coefficients (to avoid boundary effects)
Ly = (block(1)-1)/2;		% block(1) and block(2) must be odd!
Lx = (block(2)-1)/2;
if (Ly~=floor(Ly))||(Lx~=floor(Lx))
   error('Spatial dimensions of neighborhood must be odd!');
end   
nexp = nblv*nblh;	
N = prod(block)
foo = zeros(nexp,N);

% Compute covariance of noise from 'noise'
n = 0;
for ny = -Ly:Ly	% spatial neighbors
	for nx = -Lx:Lx
		n = n + 1;
		foo = shift(noise(:,:,1),[ny nx]);
		foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
		W(:,n) = vector(foo);
	end
end

C_w = innerProd(W)/nexp;

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

%Error vector
err = ones(1,3);
mean_error = mean(err);


%Initialising the covariance matrices
Cov_k = zeros(patchsize^2, patchsize^2, k);
Cov_k_n1 = zeros(patchsize^2, patchsize^2, k);
covcount = 0;

n = 0;
for ny = -Ly:Ly	% spatial neighbors
	for nx = -Lx:Lx
		n = n + 1;
		foo = shift(im(:,:,1),[ny nx]);
		foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
		W(:,n) = vector(foo);
	end
end

temp = innerProd(W)/nexp;

for i = 1:k
    Cov_k(:,:,i) = temp;
end

% for i = 1:patchsize:h-patchsize+1
%     for j = 1:patchsize:l-patchsize+1
%         h_end = i+patchsize-1;
%         l_end = j+patchsize-1;
%         patch = paddedim(i:h_end , j:l_end);
%         covtemp = cov(patch);
%         Cov_k(:, :, covcount+1) = covtemp ;
%         covcount = covcount + 1;
%         if covcount == k
%             break
%         end
%     end
%     if covcount == k
%         break
%     end
% end   


%EM Algorithm implementation
%Estimating the parameters

while( mean_error > 0.001)
    
    %Calculating other probabilities that are required
    %Calculating p(y(m)/k)
    p_ym_k = zeros(m, k);
    p_ym = zeros(1,m);
    p_ym_k_z = zeros(m,k,numz);
    for i = 1:m
        for j = 1:k
            for q = 1:numz
                p_ym_k_z(i,j,q) = gaussian_dist( patchmatrix(1:patchlen,i), zi(q)*Cov_k(:,:,j)+ C_w, numel(Cov_k(:,:,j)));
            end
            squeezed_array = p_ym_k_z(i,j,:);
            squeezed_array = squeezed_array(:,:);
            p_ym_k(i,j) = sum(times( squeezed_array, p_z));
        end
        p_ym(i) = prob_y( p_ym_k( i, :), P_k_n);
    end

    %Equation 1:
    P_k_n1 =  times ( (P_k_n).', ( sum(p_ym_k)./ p_ym )) / m;
    
    
    %Equation 2:
    intermediate = zeroes(numz, k);
    for i = 1:numz
        for j = 1:k
            intermediate(i,j) = sum( times( p_ym_k_z(:,i,j), p_ym));
        end
    end
    p_k_z_n1 = times (p_k_z_n, intermediate)./m;
    p_k_z_n1 = bsxfun(@rdivide, p_k_z_n1, P_k_n1);
     
    
    %Equation 3:
    for i= 1:k
        matsum = zeroes(patchsize, patchsize);
        for j = 1:numz
            matsum_k_z = zeroes(patchsize, patchsize);
            coeff_sum = 0;
            for k = 1:m
                y_m = reshape( patchmatrix(:,k), patchsize, patchsize);
                coeff = p_ym(k)* p_ym_k_z(k,i,j);
                matsum_k_z = matsum_k_z + coeff * ( y_m * y_m.' );
                coeff_sum = coeff_sum + coeff;
            end
            matsum = matsum + (matsum_k_z / coeff_sum) * p_z(j);
        end
        Cov_k_n1(i) = matsum - C_w; 
    end    
      
    
    %Normalise and find mean square error
    normP_k_n = P_k_n - min(P_k_n(:));
    normP_k_n = P_k_n ./ max(normP_k_n(:));
    normP_k_n1 = P_k_n1 - min(P_k_n1(:));
    normP_k_n1 = P_k_n1 ./ max(normP_k_n1(:));
    
    err(1) = immse(normP_k_n, normP_k_n1);
    
    normp_k_z_n = p_k_z_n - min(p_k_z_n(:));
    normp_k_z_n = p_k_z_n ./ max(normp_k_z_n(:));
    normp_k_z_n1= p_k_z_n1 - min(p_k_z_n1(:));
    normp_k_z_n1 = p_k_z_n1 ./ max(normp_k_z_n1(:));
    
    err(2) = immse(normp_k_z_n, normp_k_z_n1);
    
    normCov_k = Cov_k - min(Cov_k(:));
    normCov_k = Cov_k ./ max(normCov_k(:));
    normCov_k_n1 = Cov_k_n1 - min(Cov_k_n1(:));
    normCov_k_n1 = Cov_k_n1 ./ max(normCov_k_n1(:));
    
    err(3) = immse( normCov_k, normCov_k_n1);
    
    mean_error = mean(err);
    
    %Update the values
    P_k_n = P_k_n1;
    p_k_z_n = p_k_z_n1;
    Cov_k = Cov_k_n1;
    
    
end




















