function [mean_x,pdf_x]=mean_from_cdf(x,cdf_x)
pdf_x=cdf_x*0;
for i=2:length(x)
    pdf_x(i)=(cdf_x(i)-cdf_x(i -1))./(x(2)-x(1));
end
mean_x=trapz(x,x.*pdf_x');
