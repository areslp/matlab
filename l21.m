function [ W ] = l21( Q, lambda )
%L21 Summary of this function goes here
%   Detailed explanation goes here

W=zeros(size(Q,1),size(Q,2));

for i=1:size(Q,2)
    qi=Q(:,i);
    if lambda<norm(qi)
        W(:,i)=(norm(qi)-lambda)/norm(qi)*qi;
    else
        W(:,i)=0;
    end
end
end

