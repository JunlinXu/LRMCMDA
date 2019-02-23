 function out = GetseqM(X)
 %% GetseqM :transform the Seq similar matrix(SeqSM) to the form we need
 %    input  X:  the similar matrix
 %    output B:  the matrix need in the formula
b = sum(X,2)-1;
[m,n] = size(X);
B = zeros(m,n);
B = B-X;
for i = 1:m
    B(i,i) = b(i);
end
out = B;
end