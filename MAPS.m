function maps1=MAPS(maps)
%Get projection
[m,n]=size(maps);
maps1=zeros(m,m);
for i=1:m
    for j=1:m
        if(i~=j)
            for k=1:n
                if(maps(i,k)==1&&maps(j,k)==1)
                    maps1(i,j)=1;
                end
            end
        end
    end
end
end
