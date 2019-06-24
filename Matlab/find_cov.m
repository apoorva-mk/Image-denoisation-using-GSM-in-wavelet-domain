function [foo, Cov_matrix] = find_cov(y, block)
    [nv, nh, ~] = size(y);
    nblv = nv - block(1) + 1;
    nblh = nh - block(2) + 1;
    nexp = nblv * nblh ;
    Ly = (block(1) -1)/2;
    Lx = (block(2) -1)/2;
    N  = prod(block);
    
    n = 0;
    foo = zeros(nexp, N);
    for ny = -Ly:Ly
        for nx = -Lx:Lx
            n = n + 1;
            f = shift(y(:,:,1),[ny nx]);
            f = f(Ly+1:Ly+nblv, Lx+1:Lx+nblh);
            foo(:,n) = f(:);
        end
    end
    
    Cov_matrix = innerProd(foo)/nexp;
    
    