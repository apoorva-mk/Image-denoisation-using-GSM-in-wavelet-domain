clear all
clc

im0 = imread('images/barco.png');
% im0 = imresize(im0,2);
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


%Daubechies Wavelet decomposition
[cA1,cH1,cV1,cD1] = dwt2(im0,'db2');

%EM Implementation
patchsize = 5;

im = cV1;
im0 = im;
[l,h]= size(im);

%padding with zeros to make it divisible by patch size
lpad = patchsize - mod(l,patchsize);
hpad = patchsize - mod(h,patchsize);
paddedim = padarray(im, [lpad,hpad], 'post'); 

%Vectorising the neighbourhoods
%Extracting patches and vectorising the patches and making a new matrix
patchmatrix = [];
patchstep = 1;


for i = 1:patchstep:(h-patchsize+1)
    for j = 1:patchstep:(l-patchsize+1)     
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
block = [patchsize,patchsize];
nblv = nv-block(1)+1;	% Discard the outer coefficients 
nblh = nh-block(2)+1;   % for the reference (centrral) coefficients (to avoid boundary effects)
Ly = (block(1)-1)/2;		% block(1) and block(2) must be odd!
Lx = (block(2)-1)/2;
if (Ly~=floor(Ly))||(Lx~=floor(Lx))
   error('Spatial dimensions of neighborhood must be odd!');
end   
nexp = nblv*nblh;	
N = prod(block);
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
[Q,L] = eig(C_w);
% correct possible negative eigenvalues, without changing the overall variance
L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
C_w = Q*L*Q';

%Initialisation of proabability vectors 
k = 10;
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

%padding with zeros to make it divisible by patch size
lpad2 = k - mod(l,k);
hpad2 = k - mod(h,k);
paddedim2 = padarray(im, [lpad2,hpad2], 'post'); 

step = size(paddedim2)/k;



for i = 1:k
    temp = compute_cov( paddedim2((i-1)*step(1)+1:i*step ,:),patchsize)-C_w;
    [Q,L] = eig(temp);
    L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
    Cov_k(:,:,i) = Q*L*Q';
end

%EM Algorithm implementation
%Estimating the parameters
iter = 1;
while(iter<11)%( mean_error > 0.001)
    tic
    disp(['ITERATION NO. : ',num2str(iter)]);
    iter = iter+1;
        
    %Calculating other probabilities that are required
    %Calculating p(y(m)/k)
    p_ym_k = zeros(m, k);
    p_ym = zeros(1,m);
    p_ym_k_z = zeros(m,k,numz);

    for i = 1:k
        for j = 1:numz
            p_ym_k_z(:,i,j) = gaussian_dist( patchmatrix, zi(j)*Cov_k(:,:,i)+ C_w, patchlen);
        end
    end  
    
    for i = 1:numz
        p_ym_k = p_ym_k + p_ym_k_z(:,:,i)*p_z(i);
    end
    
    p_ym = p_ym_k * P_k_n';
    
    p_ym(p_ym == 0) = eps;
    
        

    %Equation 1:
    P_k_n1 =  times ( (P_k_n), ( sum(p_ym_k./ p_ym))) / m;
    
    
    %Equation 2:
    intermediate = reshape(p_ym_k_z,m, numz*k);
    intermediate = (intermediate./ p_ym);
    intermediate = sum(intermediate)./m;
    p_z_k_n_reshaped = reshape(p_z_k_n,1,numz*k);
    p_z_k_n_reshaped = times(p_z_k_n_reshaped, intermediate);
    p_z_k_n1 = reshape(p_z_k_n_reshaped, numz, k);
     
    
    %Equation 3:
    for i= 1:k
        matsum = zeros(patchlen, patchlen);
        for j = 1:numz
            matsum_k_z = zeros(patchlen, patchlen);
            coeff_sum = 0;
            for r = 1:m
                y_m = ( patchmatrix(:,r));%, patchsize, patchsize);
                coeff = p_ym_k_z(r,i,j)/p_ym(r);
                matsum_k_z = matsum_k_z + coeff * ( y_m * y_m.' );
                coeff_sum = coeff_sum + coeff;
            end
            matsum = matsum + (matsum_k_z / coeff_sum) * p_z(j);
        end
        Cov_k_n1(:,:,i) = real(matsum - C_w);
%         Cov_k_n1(:,:,i)
        [Q,L] = eig(Cov_k_n1(:,:,i));
        % correct possible negative eigenvalues, without changing the overall variance
        L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
        Cov_k_n1(:,:,i) = Q*L*Q';
    end
    
    
    %Normalise and find mean square error
    
    if iter == 2
        normP_k_n = P_k_n./max(P_k_n);
        
    else
        normP_k_n = P_k_n - min(P_k_n(:));
        normP_k_n = P_k_n ./ max(normP_k_n(:));
    end
        normP_k_n1 = P_k_n1 - min(P_k_n1(:));
        normP_k_n1 = P_k_n1 ./ max(normP_k_n1(:));
    
    err(1) = immse(normP_k_n, normP_k_n1);
    
    if iter==2
        normp_z_k_n = p_z_k_n./max(p_z_k_n);
    else
        normp_z_k_n = p_z_k_n - min(p_z_k_n(:));
        normp_z_k_n = p_z_k_n ./ max(normp_z_k_n(:));
    end
    normp_z_k_n1= p_z_k_n1 - min(p_z_k_n1(:));
    normp_z_k_n1 = p_z_k_n1 ./ max(normp_z_k_n1(:));
    
    err(2) = immse(normp_z_k_n, normp_z_k_n1);
    
    normCov_k = Cov_k - min(Cov_k(:));
    normCov_k = Cov_k ./ max(normCov_k(:));
    normCov_k_n1 = Cov_k_n1 - min(Cov_k_n1(:));
    normCov_k_n1 = Cov_k_n1 ./ max(normCov_k_n1(:));
    
    err(3) = immse( normCov_k, normCov_k_n1);
    
    mean_error = mean(err);
    
    %Calculating log-likelihood function
    loglike(iter-1) = sum(log(p_ym_k * P_k_n1'));
    
     %Update the values
     P_k_n = P_k_n1;
     p_z_k_n = p_z_k_n1;
     Cov_k = Cov_k_n1;  
     toc
end




















