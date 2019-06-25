
function [ x_hat] = denoising_utility(im,k,patchsize,noise)

%%
% initializations

block = [patchsize, patchsize];
[l,h]= size(im);
nblv = l - block(1) + 1;
nblh = h - block(2) + 1;
nexp = nblv * nblh ;
Ly = (block(1) -1)/2;
Lx = (block(2) -1)/2;
N  = prod(block);


%Extracting patches and vectorising the patches and making a new matrix
patchmatrix = [];

Y = zeros(nexp, N);
n = 0;
for ny = -Ly:Ly
    for nx = -Lx:Lx
        n = n + 1;
        f = shift(im(:,:,1),[ny nx]);
        f = f(Ly+1:Ly+nblv, Lx+1:Lx+nblh);
        Y(:,n) = vector(f);
    end
end

%m stores the number of patches
patchmatrix = Y';
[ patchlen, m ] = size(patchmatrix);

%Finding C_w
W = zeros(nexp,N);

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
            disp(strcat("Creating gaussian dist ", int2str(i), " ", int2str(j)));
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
            ppp = 1;
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

figure, plot(1:iter-1,real(loglike));

%Calculating p_k_ym which is to be returned
p_k_ym = ((p_ym_k ./ P_k_n) .* p_ym)';



% denoising starts here

% aux = zeros(m,patchlen);
% for i=1:k,
%     temp_cov = Cov_k(:,:,i);
%     first_sum = zeros(m, patchlen);
%     for j=1:numz,
%         tc = zi(j)*temp_cov;
%         itc = tc + C_w;
%         itc = (itc)\eye(size(itc));
%         res = tc*itc*patchmatrix;
%         res = res' .*(p_ym_k_z(:,:,j)*p_z_k_n(j,:)');
%         res = res ./ p_ym_k(i);
%         first_sum = first_sum + res;
%     end;
%     P_ki_ym = p_ym_k(:,i)*P_k_n(i);
%     P_ki_ym = P_ki_ym./p_ym;
%     temp_sec = first_sum.*P_ki_ym;
%     aux = aux + temp_sec';
% end;

% denoising attempt 2

% just checking denoising
N = prod(block);       % size of block length X breadth


% Stores the central coefficients
cent = floor((prod(block) + 1)/2);

% They will be used for patch matrix
% size ( M X N), M:= number of patches, N:=size of patch
[~, C_w] = find_cov(noise, block);
[S, dd] = eig(C_w);
S = S*real(sqrt(dd));   % S* S' = C_w
iS = S\eye(size(S));
sig2 = mean(diag(C_w));

denoised_res = zeros(m,1);
for t=1:k
    y = im;
    [nv, nh, nb] = size(y);

    nblv = nv - block(1) + 1;   % Discarding the outer coefficients
    nblh = nh - block(2) + 1;   %
    nexp = nblv * nblh;         % number of coefficient considered

    zM = zeros(nv, nh);    % hidden variable z
    x_hat = zeros(nv, nh); % coefficient estimation
    [Y, C_y] = find_cov(y, block);
    sy2 = mean(diag(C_y));  % observed (signal + noise) variance in
                            % the subband
                            
    C_u = Cov_k(:,:,t);
    [Q, L] = eig(iS*C_u*iS');
    la = diag(L);
    la = real(la).*(real(la)>0);
    
    V = Q'*iS*Y';
    V2 = (V.^2).';
    M = S*Q;
    m = M(cent);
    
    % Compute p(Y|log(z))
    
    % since we are usig non-informative prior
    lzmin = -20.5;
    lzmax =   3.5;
    step  =   0.5;
    
    lzi = lzmin:step:lzmax;
    nsamp_z = length(lzi);
    zi = exp(lzi);
    
    
    laz = la*zi;
    p_lz = zeros(nexp, nsamp_z);
    mu_x = zeros(nexp, nsamp_z);
    
    
    %%
    pg1_lz = 1./sqrt(prod(1 + laz,1));	% normalization term (depends on z, but not on Y)
    aux = exp(-0.5*V2*(1./(1+laz)));
    p_lz = aux*diag(pg1_lz);				% That gives us the conditional Gaussian density values
    										% for the observed samples and the considered samples of z
    % Compute mu_x(z) = E{x|log(z),Y}
    aux = diag(m)*(laz./(1 + laz));	% Remember: laz = la*zi
    mu_x = V.'*aux;	
    %%
    [foo, ind] = max(p_lz.');	% We use ML estimation of z only for the boundaries.
    clear foo
    if prod(size(ind)) == 0,
        z = ones(1,size(ind,2));
    else
        z = zi(ind).';				
    end
    % For boundary handling

    uv=1+Ly;
    lh=1+Lx;
    dv=nblv+Ly;
    rh=nblh+Lx;
    ul1=ones(uv,lh);
    u1=ones(uv-1,1);
    l1=ones(1,lh-1);
    ur1=ul1;
    dl1=ul1;
    dr1=ul1;
    d1=u1;
    r1=l1;

    zM(uv:dv,lh:rh) = reshape(z,nblv,nblh);

    % Propagation of the ML-estimated z to the boundaries

    % a) Corners
    zM(1:uv,1:lh)=zM(uv,lh)*ul1;
    zM(1:uv,rh:nh)=zM(uv,rh)*ur1;
    zM(dv:nv,1:lh)=zM(dv,lh)*dl1;
    zM(dv:nv,rh:nh)=zM(dv,rh)*dr1;
    % b) Bands
    zM(1:uv-1,lh+1:rh-1)=u1*zM(uv,lh+1:rh-1);
    zM(dv+1:nv,lh+1:rh-1)=d1*zM(dv,lh+1:rh-1);
    zM(uv+1:dv-1,1:lh-1)=zM(uv+1:dv-1,lh)*l1;
    zM(uv+1:dv-1,rh+1:nh)=zM(uv+1:dv-1,rh)*r1;
    %%

    
    p_z = ones(nsamp_z,1);    % Flat log-prior (non-informative for GSM)
    p_z = p_z/sum(p_z);

    % Compute p(log(z)|Y) from p(Y|log(z)) and p(log(z)) (Bayes Rule)

    p_lz_y = p_lz*diag(p_z);
    clear p_lz
    aux = sum(p_lz_y, 2);
    if any(aux==0),
        foo = aux==0;
        p_lz_y = repmat(~foo,1,nsamp_z).*p_lz_y./repmat(aux + foo,1,nsamp_z)...
            + repmat(foo,1,nsamp_z).*repmat(p_z',nexp,1); 	% Normalizing: p(log(z)|Y)
    else
        p_lz_y = p_lz_y./repmat(aux,1,nsamp_z); 	% Normalizing: p(log(z)|Y)
    end
    % Compute E{x|Y} = int_log(z){ E{x|log(z),Y} p(log(z)|Y) d(log(z)) }

    aux = sum(mu_x.*p_lz_y, 2);
    denoised_res = denoised_res + p_k_ym(t,:)'.*aux;

end

x_hat(1+Ly:nblv+Ly,1+Lx:nblh+Lx) = reshape(aux,nblv,nblh);

