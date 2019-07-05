function [ Cov_k, P_k_ym ] = em_(im, k, patchsize, noise)
% TODOs:
% 1. Remove the debug parts when not required
% 2. Change the ym initialization according to assumptions
% 3. Add more vectorizations to speed things up
% 4. Add the loglike results

%% Initializing step
block = [patchsize, patchsize];

[Y, ~] = find_cov(im, block);
[~, C_w] = find_cov(noise, block);
patchmatrix = Y';
m = size(patchmatrix,2);

% Change to a number comprehensible by MATLAB
% remove if things work out w/o it
% patchmatrix(patchmatrix<eps) = eps;

clear Y;

% Making C_w positive definite
[Q, L] = eig(C_w);
L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
C_w = Q*L*Q';

% Initializing z values
lzmin = -20.5;
lzmax = 3.5;
step = 2;
lzi = lzmin:step:lzmax;
numz = length(lzi);
zi = exp(lzi);
p_z = ones(1, numz)/numz;

% Initializing the cov_matrices
Cov_k = zeros(patchsize^2, patchsize^2, k);
Cov_k_n1 = zeros(patchsize^2, patchsize^2, k);

for i=1:k
    % Here we have used images after circular-shifting it
    temp = compute_cov(shift(im, [i*20, i*30]), patchsize);
    [Q, L] = eig(temp);
    L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
    Cov_k(:,:,i) = Q*L*Q';
end

%% Initializing probabilities values
P_k_n = ones(k,1)/k;
p_z_k_n = ones(numz,k)/numz;

%% EM loop

% calculation of ym * ym.T
disp("Calculating ym * ym.T");
ymymT = zeros(m, patchsize^2, patchsize^2); % dim: M X (patch)^2 X (patch)^2
ym_ = patchmatrix';
for p=1:patchsize^2
    ymymT(:,:,p) = ym_ .* ym_(:,p);
end
ymymT = reshape(ymymT, m, []); % flattened the 2nd and 3rd dims.

clear ym_;

iter = 1;
disp("E M now starts"); % DEBUG

while(iter < 11)
    tic
    disp(['Iteration No. :', num2str(iter)]); % DEBUG
    iter = iter + 1;
    
    % Initialization of probabilities
    p_ym_k = zeros(m,k);
    p_ym_k_z = zeros(m,k,numz);
    
    % p(ym | k, z) = Gaussian_dist( z * Ck + Cw, mean=0)
    disp('Calculating p(ym | k, z) using Gaussian_dist'); % DEBUG
    for i = 1:k 
        for j = 1:numz
            tt = gaussian_dist( patchmatrix, zi(j)*Cov_k(:,:,i)+ C_w, patchsize^2);
            % An attempt to remove NaN from the results upfront
            % taking gaussian dist to always return positive values
            % Change small numbers to a number comprehensible by the
            % MATLAB w/o giving NaN or Inf
            % Remove if things work out w/o it
%             tt(tt<eps) = eps; 
            
            p_ym_k_z(:,i,j) = tt;
        end
    end  
    
    % p(ym | k) = Integrate(p(ym | k, z) * pk(z), dz
    disp('Calculating p(ym | k)'); % DEBUG
    for i = 1:numz
        p_ym_k = p_ym_k + p_ym_k_z(:,:,i).*p_z_k_n(i,:);
    end
    
    % Putting an assumption that all ym are equally probable to be
    % selected so we are taking p(ym) = 1/m
    % It can be changed here if required
    p_ym = ones(m, 1) * (1/m);
    
    
    % P(k | ym) = (p(ym | k) * P(k)) / p(ym)
    P_k_ym = P_k_n .* p_ym_k';
    P_k_ym = P_k_ym ./ p_ym';
    
    % Now we have all the required probabilities and we can start
    % implementing the Equations
    
    % Equation 1:
    % Updates the P(k) for each iteration
    % P(k)' = P(k) * (1/M)
    %          * SUM((m:1,M), p(ym | k)/SUM((j:1,k), p(ym | j)*P(j)))
    
    disp('Equation 1: being done') % DEBUG
    temp = p_ym_k * P_k_n;
    temp = p_ym_k ./ temp;
    temp = sum(temp);
    P_k_n1 = (P_k_n .* temp') / m;
    
    
    % Equation 2:
    % Updates the pk(z) for each iteration
    
    disp('Equation 2: being done') % DEBUG
    P_inv_ym = P_k_ym ./ p_ym_k';  % equivalent to P(k) / p(ym)
    
    temp = zeros(size(p_ym_k_z));
    for pp=1:numz
        temp(:,:,pp) = p_ym_k_z(:,:,pp).* P_inv_ym';
    end
    
    temp = reshape(sum(temp),k,numz);
    temp = temp ./ P_k_n;
    p_z_k_n1 = (p_z_k_n .* temp')/m; % Updated here
    
    % Equation 3:
    % Updates the covariance matrices in each iteration
    disp('Equation 3: being done') % DEBUG
    for i = 1:k
        matsum = zeros(patchsize^2, patchsize^2);
        for j=1:numz
            temp = (p_ym_k_z(:,i,j)*p_z_k_n(j,i))./p_ym_k(:,i);
            temp = P_k_ym(i,:)'.*temp;
            temp1 = sum(ymymT.*temp);
            C_y_k_z = temp1/sum(temp);
            matsum = matsum + reshape(C_y_k_z, patchsize^2, patchsize^2)*p_z(j);
        end
        Cov_k_n1(:,:,i) = real(matsum - C_w);
        
        % Check for failure - Only for debugging - Remove if not req.
        % Checks for any unwanted NaN appearing in COV Matrices
        if any(any(isnan(Cov_k_n1(:,:,i))))
            disp("Unwanted NaNs appeared in Cov matrices");
        end
        
        % Make the resultant cov matrix positive definite
        [Q, L] = eig(Cov_k_n1(:,:,i));
        L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
        Cov_k_n1(:,:,i) = Q*L*Q'; 
    end
    
    % Update the values
    P_k_n = P_k_n1;
    p_z_k_n = p_z_k_n1;
    Cov_k = Cov_k_n1;
    toc
end
