function [ minmat ] = elementwisecellmin( STSnorm,Np,num_exp )
%takes the elementwise minimum of the cell array STS
minn = ones(Np,Np,2);
for i=1:num_exp%1587
    Z = STSnorm{i};
    minn(:,:,2)=Z;
    minn(:,:,1)=min(minn,[],3);
    end
    minmat = minn(:,:,1);
end

