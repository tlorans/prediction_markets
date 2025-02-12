function [d_series]=trading_order(t,n)
m=floor(t/n);
d_series=ones(t,1);
for i=1:m
    for j=1:n
        d_series((i-1)*n+j,1)=j;
    end
end
counter=0;
for i=n*m+1:t
    counter=counter+1;
    d_series(i,1)=counter;
end
