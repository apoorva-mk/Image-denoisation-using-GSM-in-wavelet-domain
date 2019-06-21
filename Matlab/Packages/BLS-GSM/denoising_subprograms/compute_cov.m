function temp = compute_cov(im0, patchsize)


[nv,nh,nb] = size(im0);
block = [patchsize, patchsize];
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


n = 0;
for ny = -Ly:Ly	% spatial neighbors
	for nx = -Lx:Lx
		n = n + 1;
		foo = shift(im0(:,:,1),[ny nx]);
		foo = foo(Ly+1:Ly+nblv,Lx+1:Lx+nblh);
		W(:,n) = vectorise(foo);
	end
end

temp = innerProd(W)/nexp;

[Q,L] = eig(temp);
% correct possible negative eigenvalues, without changing the overall variance
L = diag(diag(L).*(diag(L)>0))*sum(diag(L))/(sum(diag(L).*(diag(L)>0))+(sum(diag(L).*(diag(L)>0))==0));
temp = Q*L*Q';
