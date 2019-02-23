function out = getoptT(X,W,Y,Z,S,M_E,E,m0,lam1,lam2,diseasesimilarities,lam4,mfs)
norm2WZ = norm(W,'fro')^2 + norm(Z,'fro')^2;
f(1) = F_t(X, Y,S,M_E,E,m0,lam1,lam2,diseasesimilarities,lam4,mfs) ;

t = -1e-1 ;
for i = 1:20
        f(i+1) = F_t(X+t*W,Y+t*Z,S,M_E,E,m0,lam1,lam2,diseasesimilarities,lam4,mfs) ;

        if( f(i+1) - f(1) <= .5*(t)*norm2WZ )
            out = t ;
            return;
        end
        t = t/2 ;
end
out = t ;
end