function out = getoptS1(X,S,Y,M_E,E,lam2,lam4,M,D)
M= GetseqM(M);
D= GetseqM(D);
for step=1:5
 S1=S;
YY=Y'*Y;
XX=X'*X;
A=X'*((X*S*Y').*E)*Y;
B=X'*M*X*S*YY;
BB=XX*S*Y'*D'*Y;
C=X'*(M_E.*E)*Y;
S = S.*(C./(A+2*lam2*B+2*lam4*BB)) ;
 error =mean(mean(abs(S-S1)))/mean(mean(S1));
if error< 10^(-7)
            break;
        end       
end
out =S;
end