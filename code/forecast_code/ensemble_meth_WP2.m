function []=ensemble_meth_WP3(Daily_in,Daily_out)
%%
useful_start=13;
useful_end=36;
WP_cap=107;
which_WP=2;
num_days=length(Daily_out);
num_points=100;
delta=length(Daily_out{1,1}.mu)-length(Daily_out{3,1}.mu);
wpp=linspace(0,1,num_points);
cdf_1=zeros(24,num_days,num_points);
cdf_2=zeros(24,num_days,num_points);
cdf_3=zeros(24,num_days,num_points);
cdf_4=zeros(24,num_days,num_points);
cdf_simple=zeros(24,num_days,num_points);
pdf_1=zeros(24,num_days,num_points);
pdf_2=zeros(24,num_days,num_points);
pdf_3=zeros(24,num_days,num_points);
pdf_4=zeros(24,num_days,num_points);
pdf_simple=zeros(24,num_days,num_points);
pdf_combined=zeros(24,num_days,num_points);
cdf_combined=zeros(24,num_days,num_points);
bid_simple=zeros(24,num_days);
W1_t_d=zeros(24,num_days);
W2_t_d=zeros(24,num_days);
W3_t_d=zeros(24,num_days);

for n=1:num_days
    for h=1:24
%cdf_1(h,n,:) = cdf('Normal',wpp,Daily_out{1,n}.mu(h+delta),Daily_out{1,n}.sigm(h+delta));
cdf_1(h,n,:) = cdf(Daily_out{1,n}.dis(h+delta),wpp);
pdf_1(h,n,:) = pdf(Daily_out{1,n}.dis(h+delta),wpp);
mu_1(h,n)=1-trapz(wpp,squeeze(cdf_1(h,n,:)));
cdf_2(h,n,:) = cdf(Daily_out{2,n}.dis(h+delta),wpp);
pdf_2(h,n,:) = pdf(Daily_out{2,n}.dis(h+delta),wpp);
mu_2(h,n)=1-trapz(wpp,squeeze(cdf_2(h,n,:)));
cdf_3(h,n,:) = cdf('beta',wpp,Daily_out{3,n}.a(h),Daily_out{3,n}.b(h));
pdf_3(h,n,:) = pdf('beta',wpp,Daily_out{3,n}.a(h),Daily_out{3,n}.b(h));
mu_3(h,n)=1-trapz(wpp,squeeze(cdf_3(h,n,:)));
cdf_4(h,n,:)=Daily_out{4,n}.cdf(h,:); 
mu_4(h,n)=1-trapz(wpp,squeeze(cdf_4(h,n,:)));
cdf_simple(h,n,:)=(cdf_1(h,n,:)+cdf_2(h,n,:)+cdf_3(h,n,:))/3;
mu_5(h,n)=1-trapz(wpp,squeeze(cdf_simple(h,n,:)));
real_wp(h,n)=Daily_in{4,n}(h+delta);
%*** bidding based on ***
which_val=Daily_in{7,n}(h);
distance_simple=(cdf_simple(h,n,:)-which_val).^2;
[~,find_val_simple]=min(distance_simple);
bid_simple(h,n)=wpp(find_val_simple);
    end
[index_error_ens(:,n),dev_ens(:,n)]=assessment2(WP_cap,Daily_in{4,n},bid_simple(:,n),Daily_in{6,n},useful_start,useful_end);
%a_1=Daily_out{1,n}.mu;
%b_1=Daily_out{1,n}.sigm;
%a_2=Daily_out{2,n}.mu;
%b_2=Daily_out{2,n}.sigm;
a_3=Daily_out{3,n}.a;
b_3=Daily_out{3,n}.b;
pd1=Daily_out{1,n}.dis(h+delta);
pd2=Daily_out{2,n}.dis(h+delta);
first_h_day=(n+hist_step/24)*24-(useful_start-1)-hist_step;
%[W_t]=wights_combined_temp(a_1',b_1',a_2,b_2,a_3,b_3,which_WP,first_h_day,WP,hist_step);
[W_t]=wights_combined_temp(pd1,pd2,a_3,b_3,which_WP,first_h_day,WP,hist_step);
W1_t_d(:,n)=W_t(:,1);
W2_t_d(:,n)=W_t(:,2);
W3_t_d(:,n)=W_t(:,3);
end
%%
for n=1:num_days
    for h=1:24
weight1(h,n)=W1_t_d(h,n)/(W1_t_d(h,n)+W2_t_d(h,n)+W3_t_d(h,n));
weight2(h,n)=W2_t_d(h,n)/(W1_t_d(h,n)+W2_t_d(h,n)+W3_t_d(h,n));
weight3(h,n)=W3_t_d(h,n)/(W1_t_d(h,n)+W2_t_d(h,n)+W3_t_d(h,n));
mu_6(h,n)=weight1(h,n)*mu_1(h,n)+weight2(h,n)*mu_2(h,n)+weight3(h,n)*mu_3(h,n);
pdf_combined(h,n,:)=weight1(h,n)*pdf_1(h,n,:)+weight2(h,n)*pdf_2(h,n,:)+weight3(h,n)*pdf_3(h,n,:);
cdf_combined(h,n,:)=make_cdf(wpp,pdf_combined(h,n,:));
which_val=Daily_in{7,n}(h);
distance_combined=(cdf_combined(h,n,:)-which_val).^2;
[~,find_val_combined]=min(distance_combined);
bid_combined(h,n)=wpp(find_val_combined);
    end
    [index_error_combined(:,n),dev_combined(:,n)]=assessment2(WP_cap,Daily_in{4,n},bid_combined(:,n),Daily_in{6,n},useful_start,useful_end);
end
%%

% figure()
% plot(squeeze(pdf(Daily_out{1,14}.dis(1+delta),wpp)))
% hold on
% plot(squeeze(pdf(Daily_out{2,14}.dis(1+delta),wpp)))
% plot(squeeze(betapdf(wpp,Daily_out{3,14}.a(1),Daily_out{3,14}.b(1))))
% plot(squeeze(pdf_combined(1,14,:)))
% legend
% 
% 
% 
% 
% %%
% figure()
% plot(weight1(:,14))
% hold on
% plot(weight2(:,14))
% plot(weight3(:,14))
% legend
%%
Daily_RMSE=100*sqrt(sum((daily_deviations.^2)))/sqrt(24);
%********************
%%
total_deviation2=three_to_two(daily_deviations);
total_dev_ens=reshape(dev_ens,[(24*num_days),1]);
total_dev_combined=reshape(dev_combined,[(24*num_days),1]);
sum_dev=[sum(total_deviation2),sum(total_dev_ens),sum(total_dev_combined)];
sum_abs_dev=[sum(abs(total_deviation2)),sum(abs(total_dev_ens)),sum(abs(total_dev_combined))];
total_MBE=[100*sum(total_deviation2)/(length(days)*24),100*sum(total_dev_ens)/(length(days)*24),100*sum(total_dev_combined)/(length(days)*24)];
total_MAE=[100*sum(abs(total_deviation2))/(length(days)*24),100*sum(abs(total_dev_ens))/(length(days)*24),100*sum(abs(total_dev_combined))/(length(days)*24)]
total_MSE=[100*sum((total_deviation2.^2))/(length(days)*24),100*sum((total_dev_ens.^2))/(length(days)*24),100*sum((total_dev_combined.^2))/(length(days)*24)];
total_RMSE2=[100*sqrt(sum((total_deviation2.^2))/(length(days)*24)),100*sqrt(sum((total_dev_ens.^2))/(length(days)*24)),100*sqrt(sum((total_dev_combined.^2))/(length(days)*24))];
total_error2=[sum(daily_error2'),sum(index_error_ens(2,:)),sum(index_error_combined(2,:))];
total_error3=[sum(daily_error3'),sum(index_error_ens(3,:)),sum(index_error_combined(3,:))];
%%
%*******ensemble mean forecasting
total_mu1=reshape(mu_1,[(24*num_days),1]);
total_mu2=reshape(mu_2,[(24*num_days),1]);
total_mu3=reshape(mu_3,[(24*num_days),1]);
total_mu4=reshape(mu_4,[(24*num_days),1]);
total_mu5=reshape(mu_5,[(24*num_days),1]);
total_mu6=reshape(mu_6,[(24*num_days),1]);
total_real=reshape(real_wp,[(24*num_days),1]);
Devs_mu=total_real-[total_mu1,total_mu2,total_mu3,total_mu4,total_mu5,total_mu6];
ensemble_MAE=100*sum(abs(Devs_mu))/(length(days)*24);
ensemble_RMSE=100*sqrt(sum((Devs_mu.^2))/(length(days)*24));
%%
%*********Daily Indexes******
daily_dev1=[real_wp-mu_1];
Daily_RMSE_dev1=100*sqrt(sum((daily_dev1.^2)))/sqrt(24);

daily_dev2=[real_wp-mu_2];
Daily_RMSE_dev2=100*sqrt(sum((daily_dev2.^2)))/sqrt(24);

daily_dev3=[real_wp-mu_3];
Daily_RMSE_dev3=100*sqrt(sum((daily_dev3.^2)))/sqrt(24);

daily_dev4=[real_wp-mu_4];
Daily_RMSE_dev4=100*sqrt(sum((daily_dev4.^2)))/sqrt(24);

daily_dev5=[real_wp-mu_5];
Daily_RMSE_dev5=100*sqrt(sum((daily_dev5.^2)))/sqrt(24);

RMSE_daily_matrice=[Daily_RMSE_dev1;Daily_RMSE_dev2;Daily_RMSE_dev3;Daily_RMSE_dev4;Daily_RMSE_dev5];
%%
% cnt=0;
% for i=1:num_days
%     if min(RMSE_daily_matrice(:,i))==RMSE_daily_matrice(3,i)
%         cnt=cnt+1;
%     end
% end
% 

%%
% close all
% figure()
% plot([1:24],bid_simple(:,75))
% hold on
% plot([1:24],Daily_out{4,75}.bid)
% legend
% %%
% %close all
% figure()
% plot(wpp,squeeze(cdf_1(15,75,:)))
% hold on
% plot(wpp,squeeze(cdf_2(15,75,:)))
% plot(wpp,squeeze(cdf_3(15,75,:)))
% plot(wpp,squeeze(cdf_4(15,75,:)))
% plot(wpp,squeeze(cdf_simple(15,75,:)))
% legend



