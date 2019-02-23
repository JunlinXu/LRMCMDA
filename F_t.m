function out = F_t(X,Y,S,M_E,E,m0,lam1, lam2,diseasesimilarities,lam4,mfs)
[n,r] = size(X) ;
[m,r] = size(Y) ;
s1 = diseasesimilarities;
s=mfs;
a = X*S*Y';
[m,n] = size(s);
[m1,n1] = size(s1);
a1=Y*S*X';
Sum = 0;
Sum1=0;
for i = 2:m
     for j = 1:i-1
         Sum = Sum+s(i,j)*norm((a(i,:)-a(j,:)),2).^2;
     end
end
for i = 2:n1
     for j = 1:i-1
         Sum1 = Sum1+s1(i,j)*norm((a1(:,i)-a1(:,j)),2).^2;
     end
end
out4 =lam2*Sum;
out5 =lam4*Sum1;
out1 = sum( sum( ( (X*S*Y' - M_E).*E ).^2 ) )/2 ;
out2 =  lam1*G(Y,m0,r) ;
out3 =  lam1*G(X,m0,r) ;
out = out1+out2+out3+out4+out5;
end
function out = G(X,m0,r)
z = sum(X.^2,2)/(2*m0*r) ;
y = exp( (z-1).^2 ) - 1 ;
y( find(z < 1) ) = 0 ;
out = sum(y) ;
end