%To find the probability for a given covariance matrix and zero mean)

function prob_density = gaussian_dist( x, covariance_matrix, neighbourhood_size)

num = exp (x.' * covariance_matrix * x * (-0.5)) ;
denom =  ( (2*(pi))^(neighbourhood_size/2) * (det(covariance_matrix)^(0.5)) );
prob_density = diag(num)/denom;


