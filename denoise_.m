function [denoised_y] = denoise_(y, block, noise, cov_matrices, prob_dist)

[nv, nh, ~] = size(y);
nblv = nv - block(1) + 1;
nblh = nh - block(2) + 1;
nexp = nblv * nblh;
Ly = (block(1) -1)/2;
Lx = (block(2) -1)/2;

N = prod(block);
cent = floor((prod(block) + 1)/2);

[~, C_w] = find_cov(noise, block);
[S, dd] = eig(C_w);
S = S*real(sqrt(dd));
iS = S\eye(size(S));
sig2 = mean(diag(C_w));

[Y, C_y] = find_cov(y, block);
sy2 = mean(diag(C_y));


m = nexp;


k = size(cov_matrices, 3);

denoised_y = zeros(size(y));


denoised_res = zeros(m,1);

for t=1:k
    zM = zeros(nv, nh);
    x_hat = zeros(nv, nh);
    
    C_u = cov_matrices(:,:,t);
    sx2 = sy2 - sig2;
    sx2 = sx2.*(sx2>0);
    [Q, L] = eig(iS*C_u*iS');
    
    la = diag(L);
    la = real(la).*(real(la)>0);
    
    V = Q'*iS*Y';
    V2 = (V.^2).';
    M = S*Q;
    m = M(cent,:);
    
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
    
    x_hat = y(:,:,1).*(sx2*zM)./(sx2*zM + sig2);
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
    denoised_res = denoised_res + prob_dist(t,:)'.*aux;
    
end

denoised_y(1+Ly:nblv+Ly,1+Lx:nblh+Lx) = reshape(denoised_res,nblv,nblh);

