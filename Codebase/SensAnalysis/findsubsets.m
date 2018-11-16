function [ subsets ] = findsubsets( A,params )
%returns a vector where the fist column's parameters are a subset of the
%corresponding parameter in the second columns
count=0;
for i=1:22
    for j=1:22
        if j~=i
        if all(ismember(A(i,:),A(j,:)))
            %disp(params(i))
            %disp('is a subset of')
            count=count+1;
            %disp(params(j))
        end
        end
    end
end
subsets=zeros(count,2);
for i=1:22
    for j=1:22
        if j~=i
        if all(ismember(A(i,:),A(j,:)))
            %disp('is a subset of')
            count=count+1;
            subsets(count,1)=(i)
            subsets(count,2)=(j)
        end
        end
    end
end
end
