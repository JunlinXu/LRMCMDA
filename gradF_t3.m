function [W, Z] = gradF_t3(X,Y,S,M_E,E,m0,lam1,lam2,lam4,dss1,mfs)
[n, r] = size(X);
[m, r] = size(Y);
XS = X*S;
YS = Y*S';
 XSY = XS*Y';
 A= GetseqM(dss1);
B= GetseqM(mfs);
G3y =YS*X'*B*XS;
G3x = B*XSY*YS;
G4x =XSY*A*YS;
G4y =A*YS*X'*XS;
Qx = X'* ( (M_E - XSY).*E )*YS /n;
Qy = Y'* ( (M_E - XSY).*E )'*XS /m;
W = ( (XSY - M_E).*E )*YS+ lam1*Gp(X,m0,r)+4*lam2*G3x+4*lam4*G4x;
Z = ( (XSY - M_E).*E )'*XS+ lam1*Gp(Y,m0,r)+4*lam2*G3y+4*lam4*G4y;
end

function out = Gp(X,m0,r)
z = sum(X.^2,2) /(2*m0*r) ;
z = 2*exp( (z-1).^2 ).*(z-1) ;
z( find(z<0) ) = 0;
out = X.*repmat(z,1,r) / (m0*r) ;
end



