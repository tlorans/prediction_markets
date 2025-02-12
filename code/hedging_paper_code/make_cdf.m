function [x_range,y,z]=make_cdf(x)
min_x=min(x);
max_x=max(x);
x_range=linspace(min_x,max_x,1e3)';
y=zeros(size(x_range,1),1);
for i=1:size(x_range)-1
    y(i+1)=sum( (x>=x_range(i))&(x<x_range(i+1)) );
end
sum(y);
y=y/sum(y);
z=zeros(size(y,1),1);
for i=2:size(y,1)
    z(i)=z(i-1)+y(i);
end