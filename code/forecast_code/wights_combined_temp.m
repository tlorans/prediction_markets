%function [W_t]=wights_combined_temp(a_1,b_1,a_2,b_2,a_3,b_3,which_WP,first_h_day,WP,hist_step)
function [W_t]=wights_combined_temp(pd1,pd2,a_3,b_3,which_WP,first_h_day,WP,hist_step)
%%
a1=first_h_day;
a2=first_h_day+(hist_step/3)*2;
WP_hourly=WP(which_WP).hourly_output;
WP_cap=WP(which_WP).capacity;


for t=0:35
    i=a1+t;
    k=1;
    while i<=a2   
    y_test_ML1(k,t+1)=WP_hourly(i)/WP_cap;       
    i=i+24;
    k=k+1;
    end
end
[r c]=size(y_test_ML1);
for i=1:r
    for j=1:c
        if y_test_ML1(i,j)<=0
            y_test_ML1(i,j)=0.001;
        end
    end
end

y_test_ML2=y_test_ML1(:,12:35);
W_t=zeros(24,3);
%%
for tt=1:24    
    %y_curtailed=y_test_ML2(y_test_ML2(:,tt)>0.05,tt);    
    y_curtailed=y_test_ML2(:,tt);
    x=y_curtailed;
    %f_combined=@(x,w1,w2)((w1*normpdf(x,a_1(tt+12),b_1(tt+12))+(1-w1-w2)*normpdf(x,a_2(tt+12),b_2(tt+12))+w2*betapdf(x,a_3(tt),b_3(tt))));
    %f_combined=@(x,w1,w2,w3)((w1/(w1+w2+w3))*normpdf(x,a_1(tt+12),b_1(tt+12))+(w2/(w1+w2+w3))*normpdf(x,a_2(tt+12),b_2(tt+12))+(w3/(w1+w2+w3))*betapdf(x,a_3(tt),b_3(tt)));
    opt= statset('MaxIter',1e5,'MaxFunEvals',1e5);
    f_combined=@(x,w1,w2,w3)((w1/(w1+w2+w3))*pdf(pd1,x)+(w2/(w1+w2+w3))*pdf(pd2,x)+(w3/(w1+w2+w3))*betapdf(x,a_3(tt),b_3(tt)));
    %phat=mle(x,'pdf',@(x,w1,w2)f_combined(x,w1,w2),'start',[1/3,1/3],'Options',opt,'Lowerbound',[0,0],'Upperbound',[1,1],'Optimfun','fmincon');
    phat=mle(x,'pdf',@(x,w1,w2,w3)f_combined(x,w1,w2,w3),'start',[1/3,1/3,1/3],'Options',opt,'Lowerbound',[0,0,0],'Optimfun','fmincon');
    d=phat;
    W_t(tt,1)=d(1);
    W_t(tt,2)=d(2);
    W_t(tt,3)=d(3);
end


