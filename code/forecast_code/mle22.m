function [a_beta,b_beta]=mle22(which_WP,first_h_day,WP,hist_step)
% a1=first_h_day+hist_step-forecast_step-720+48;
% a2=first_h_day+hist_step-forecast_step+48;
% a1=first_h_day+hist_step-720+48;
% a2=first_h_day+hist_step+48; 
a1=first_h_day;
a2=first_h_day+hist_step;


WP_hourly=WP(which_WP).hourly_output;
WP_cap=WP(which_WP).capacity;

%%
for t=0:35
    i=a1+t;
    k=1;
    while i<=a2   
    y_test_ML1(k,t+1)=WP_hourly(i)/WP_cap;       
    i=i+24;
    k=k+1;
    end
end
[r,c]=size(y_test_ML1);
for i=1:r
    for j=1:c
        if y_test_ML1(i,j)<=0
            y_test_ML1(i,j)=0.05;
        end
    end
end
y_test_ML2=y_test_ML1(:,12:35);
for tt=1:24
    y_curtailed=y_test_ML2(:,tt);
    phat= betafit(y_curtailed);
    d=phat;
    a_beta(tt)=d(1);
    b_beta(tt)=d(2);
    mu_beta(tt)=a_beta(tt)/(a_beta(tt)+b_beta(tt));
    beta_sigma(tt)=(a_beta(tt)*b_beta(tt))/((a_beta(tt)+b_beta(tt)+1)*(a_beta(tt)+b_beta(tt))*(a_beta(tt)+b_beta(tt)));
end
