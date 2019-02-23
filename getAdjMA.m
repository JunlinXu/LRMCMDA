function I1=getAdjMA(hou,ji)
%Construction of negative sample adjacency matrix
[m,n]=size(hou);
matDT1=zeros(m,n);
for i=1:m
    for j=1:m
        if(ji(i,j)==1&&i~=j)
            for k=1:n
                if(hou(j,k)==1&&hou(i,k)==0)
                    matDT1(i,k)=1;
                    
                end
            end
        end
        
    end
end
I=zeros(m,n);
for i=1:m
    for j=1:n
        if(matDT1(i,j)==0)
            I(i,j)=1;
        end
    end
end
I1=zeros(m,n);
for i=1:m
    for j=1:n
        if(I(i,j)==1&&hou(i,j)==0)
            I1(i,j)=1e-30;
        end
        if(I(i,j)==1&&hou(i,j)==1)
            I1(i,j)=1;
        end
    end
end
end