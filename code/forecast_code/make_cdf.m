function [cdf_x]=make_cdf(x,pdf_x)
cdf_x(1)=0;
for i=2:length(x)
    cdf_x(i)=trapz(x(1:i),pdf_x(1:i));
end