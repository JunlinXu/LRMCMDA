function [X S Y dist] = OptSpaceII(M_E,dss1,mfs,r,niter,tol,lam1,lam2,lam4)
% An algorithm for Matrix Reconstruction from a partially revealed set. 
    E=spones(M_E);
	M_E = sparse(M_E);
	[n m] = size(M_E);
    m0 = size(M_E,1);
	eps = nnz(E)/sqrt(m*n) ;
     EOrigin = E;
    rescal_param = sqrt( nnz(E) * r / norm(M_E,'fro')^2 ) ;
    M_E = M_E * rescal_param ;
fprintf(1,'Sparse SVD ...\n');
% Sparse SVD
[X0 S0 Y0] = svds(M_E,r) ;
X0 = X0*sqrt(n) ; Y0 = Y0*sqrt(m) ;
S0 = S0 /sqrt(m*n) ;
fprintf(1,'Iteration\tFit Error\n');
% Gradient Descent
X = X0;Y=Y0;S=S0;
S= getoptS1(X,S,Y,M_E,E,lam2,lam4,mfs,dss1);
matrix1=X*S*Y';
dist(1) = norm( (M_E - matrix1).*E ,'fro')/sqrt(nnz(E) )  ;
for i = 1:niter
     recMatrix = X*S*Y';
% Compute the Gradient 
	[W Z] = gradF_t3(X,Y,S,M_E,E,m0,lam1,lam2,lam4,dss1,mfs);
% Line search for the optimum jump length	
	t = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,dss1,lam4,mfs);
	X = X + t*W;
    Y = Y + t*Z;
	S = getoptS1(X,S,Y,M_E,E,lam2,lam4,mfs,dss1);
% Compute the distortion	
	dist(i+1) = norm( (recMatrix - X*S*Y').*E,'fro' )/sqrt(nnz(EOrigin));
	if( dist(i+1) < tol )
		break ;
	end
end

S = S /rescal_param ;
end