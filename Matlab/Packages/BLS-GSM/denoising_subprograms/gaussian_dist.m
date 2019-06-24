%To find the probability for a given covariance matrix and zero mean)

function prob_density = gaussian_dist( x, covariance_matrix, neighbourhood_size)
inv_cov = covariance_matrix\eye(size(covariance_matrix));
%inv_cov = covariance_matrix;
num = exp ((x' * inv_cov).* x' * (-0.5)) ;
num = sum(num,2);
denom =  ( (2*(pi))^(neighbourhood_size/2) * (det(covariance_matrix)^(0.5)) );
prob_density = num/denom;


    