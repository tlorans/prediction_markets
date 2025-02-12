clc
clear all
close all
%%
%******hourly day-ahead market price******%
%filename4 = 'da_before.xlsx';
filename5 = 'da_after.xlsx';
%da_half_before= xlsread(filename4);
da_half_after= xlsread(filename5);

%da_hourly_before=zeros(24,1);
da_hourly_after=zeros(24,1);


% for hn_bf=1:24
%     da_hourly_before(hn_bf,1)=mean(da_half_before((hn_bf-1)*2+1:hn_bf*2));
% end

for hp_af=1:24
    da_hourly_after(hp_af,1)=mean(da_half_after((hp_af-1)*2+1:hp_af*2));
end

%%
from_excel=xlsread('ImbalanceNrvPrices_201808.xls');
% MIP=from_excel((3:1442),4);
% MDP=from_excel((3:1442),5);
MIP=from_excel((1000:1442),4);
MDP=from_excel((1000:1442),5);
MDP(MDP<0)=0;
count=0;
for hd=1:length(MIP)/4
    hourly_MIP(hd,1)=mean(MIP((hd-1)*4+1:hd*4));
    hourly_MDP(hd,1)=mean(MDP((hd-1)*4+1:hd*4));
end
%******** Historical mena for forecasting **********
for ip=1:24
    daily_MID(ip,1)=mean(hourly_MIP(ip:24:(end-24+ip)));
    daily_MDP(ip,1)=mean(hourly_MDP(ip:24:(end-24+ip)));
end

%%
% filename1 = 'pos_dev.xlsx';
% filename2 = 'neg_dev.xlsx';
% lamda_pos_15min= xlsread(filename1);
% lamda_neg_15min= xlsread(filename2);
% 
% lamda_pos=zeros(24,1);
% lamda_neg=zeros(24,1);
% 
% for hp=1:24
%     lamda_pos(hp,1)=mean(lamda_pos_15min((hp-1)*4+1:hp*4));
% end
% for hn=1:24
%     lamda_neg(hn,1)=mean(lamda_neg_15min((hn-1)*4+1:hn*4));
% end
%***********
%%
lamda_pos_after=daily_MID;
lamda_neg_after=daily_MDP;
typ=lamda_pos_after/max(lamda_pos_after);
tyn=lamda_neg_after/max(lamda_neg_after);
lamda_pos_after=(1-typ).*da_hourly_after;
lamda_neg_after=(1+tyn).*da_hourly_after;
lamda_pos=da_hourly_after-lamda_pos_after;
lamda_neg=lamda_neg_after-da_hourly_after;
ratio=lamda_pos./(lamda_pos+lamda_neg);

%%
%***plot Fig .3 : electricity prices*********
color_da=[0.82 0.41 0.12];
color_pos=[1 0.41 0.71];
color_neg=[0.54 0.17 0.89];
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
stairs(da_hourly_after,'color',color_da,'linewidth',1.5)
hold on
stairs(lamda_pos,'color',color_pos,'linewidth',1.5)
stairs(lamda_neg,'color',color_neg,'linewidth',1.5)
ylabel('Electrcity Price (Euro/MWh)', 'rotation',0)
xlabel('Hour of the day')
legend ('day-ahead market clearing price','penalty price for positive deviation','penalty price for negative deviation')
box off
legend box off


%%
bid_agent1=zeros(24,1);
bid_agent2=zeros(24,1);
bid_agent3=zeros(24,1);


[y_pred_arma2,y_sd_arma2]=Arma2();
%%
%%test gpr
[y_pred2,y_sd2]=GPR2();
%%
[a_beta,b_beta]=mle22();
num_points=100;
wpp=linspace(0,1,num_points);
wind_power=0.1:0.1:0.9;
b_1=20;
p_settlement=zeros(length(wind_power),24);
for hr=1:24
cdf_Arma1 = cdf('Normal',wpp,y_pred_arma2(hr),y_sd_arma2(hr));
which_val=ratio(hr);
distance_arma=(cdf_Arma1-which_val).^2;
[~,find_val_arma]=min(distance_arma);
bid_agent1(hr)=wpp(find_val_arma);
GPR1_pd(hr) = makedist('Normal','mu',y_pred2(hr),'sigma',y_sd2(hr));
GPR1_pd_t(hr)= truncate(GPR1_pd(hr),0,1)
mu_GPR1(hr)=mean(GPR1_pd_t(hr));
sigma_GPR1(hr)=std(GPR1_pd_t(hr));
%cdf_GPR1= cdf('Normal',wpp,y_pred2(hr),y_sd2(hr));
%cdf_GPR1= cdf('Normal',wpp,mu_GPR1(hr),sigma_GPR1(hr));
cdf_GPR1= cdf(GPR1_pd_t(hr),wpp);
distance_GPR=(cdf_GPR1-which_val).^2;
[~,find_val_GPR]=min(distance_GPR);
bid_agent2(hr)=wpp(find_val_GPR);
cdf_mle1= cdf('beta',wpp,a_beta(hr),b_beta(hr));
distance_mle=(cdf_mle1-which_val).^2;
[~,find_val_mle]=min(distance_mle);
bid_agent3(hr)=wpp(find_val_mle);
% Aggregate cdf
for kb=1:length(wind_power)
which_vall=wind_power(kb);
distance=(wpp-which_vall).^2;
[~,find_vall]=min(distance);
agent_1(kb)=cdf_Arma1(find_vall);
agent_2(kb)=cdf_GPR1(find_vall);
agent_3(kb)=cdf_mle1(find_vall);
[p_settlement(kb,hr),price2(kb,hr,:)]=traders_models(agent_1(kb),agent_2(kb),agent_3(kb),b_1);
p_settled(kb)=p_settlement(kb,hr);
dfx(kb)=(exp(0/b_1).*(1/p_settled(kb) - 1))./(b_1.*(exp(0./b_1).*(1/p_settled(kb) - 1) + 1).^2);
end
cdf_aggr(hr,:)=aggregate_result(p_settled,dfx);
distance_aggr=(cdf_aggr(hr,:)-which_val).^2;
[~,find_val_aggr]=min(distance_aggr);
bid_aggr(hr)=wpp(find_val_aggr);
end
%%    
%*********plot cdfs************%
%********hour=6****************%
% figure()
% title ('hour 6')
% subplot(2,1,1)
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot(wpp,cdf('Normal',wpp,y_pred_arma2(6),y_sd_arma2(6)))
% hold on 
% plot(wpp,cdf(GPR1_pd_t(6),wpp))
% plot(wpp,cdf('beta',wpp,a_beta(6),b_beta(6)))
% plot(wpp,cdf_aggr(6,:),'k--')
% plot(wind_power,p_settlement(:,6),'O')
% legend('ARMA','GPR','MLE','aggregated','Prediction market settled price')
% xlabel('Wind power (p.u.)')
% ylabel('Cumulative Probability','rotation',0)
% legend box off
% box off
% hold off
% subplot(2,1,2)
% title ('hour 6')
% ******plot price convergence in prediction market*****%
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot([0:17],squeeze(price2(5,6,:)))
% hold on
% plot([0:17],squeeze(price2(6,6,:)))
% plot([0:17],squeeze(price2(7,6,:)))
% xlabel('Trading Round')
% ylabel(' Price in prediction market','rotation',0)
% legend('P_w < 0.5 p.u.','P_w < 0.6 p.u.','P_w < 0.7 p.u.')
% legend box off
% box off

%********hour=23****************%
% title ('hour 23')
% figure()
% subplot(2,1,1)
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot(wpp,cdf('Normal',wpp,y_pred_arma2(23),y_sd_arma2(23)))
% hold on 
% plot(wpp,cdf(GPR1_pd_t(23),wpp))
% plot(wpp,cdf('beta',wpp,a_beta(23),b_beta(23)))
% plot(wpp,cdf_aggr(23,:),'k--')
% plot(wind_power,p_settlement(:,23),'O')
% legend('ARMA','GPR','MLE','aggregated','Prediction market settled price')
% xlabel('Wind power (p.u.)')
% ylabel('Cumulative Probability','rotation',0)
% legend box off
% box off
% hold off
% subplot(2,1,2)
% %******plot price convergence in prediction market*****%
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot([0:17],squeeze(price2(5,23,:)))
% hold on
% plot([0:17],squeeze(price2(6,23,:)))
% plot([0:17],squeeze(price2(7,23,:)))
% xlabel('Trading Round')
% ylabel(' Price in prediction market','rotation',0)
% legend('P_w < 0.5 p.u.','P_w < 0.6 p.u.','P_w < 0.7 p.u.')
% legend box off
% box off
% %********hour=10****************%
%%
color_Arma=[0.78 0.08 0.52];
color_GPR=[0.12 0.56 1];
color_MLE=[0.13 0.55 0.13];
color_aggr=[0.1,0.1,0.44];
color_price=[1 0 0];
% figure()
% title ('hour 10')
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot(wpp,cdf('Normal',wpp,y_pred_arma2(10),y_sd_arma2(10)),'color',color_Arma,'linewidth',1.5)
% hold on 
% plot(wpp,cdf(GPR1_pd_t(10),wpp),'color',color_GPR,'linewidth',1.5)
% plot(wpp,cdf('beta',wpp,a_beta(10),b_beta(10)),'color',color_MLE,'linewidth',1.5)
% plot(wpp,cdf_aggr(10,:),'color',color_aggr,'linewidth',1.5,'linestyle','--')
% plot(wind_power,p_settlement(:,10),'O','color',color_price)
% legend('Agent1-ARMA','Agent2-GPR','Agent3-MLE','Aggregated in prediction market','Settled price in prediction market')
% xlabel('Wind power (p.u.)')
% ylabel('Cumulative Probability','rotation',0)
% legend box off
% box off
% hold off
% figure()
% title ('hour 10')
% %******plot price convergence in prediction market*****%
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot([0:17],squeeze(price2(5,10,:)))
% hold on
% plot([0:17],squeeze(price2(6,10,:)))
% plot([0:17],squeeze(price2(7,10,:)))
% xlabel('Trading Round')
% ylabel(' Price in prediction market','rotation',0)
% legend('If P_w < 0.5 p.u. ?','If P_w < 0.6 p.u. ?','If P_w < 0.7 p.u. ?')
% legend box off
% box off
%%
%****************next time**********************
%********hour=19****************%
figure()
title ('hour 19')
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot(wpp,cdf('Normal',wpp,y_pred_arma2(19),y_sd_arma2(19)),'color',color_Arma,'linewidth',1.5)
hold on 
plot(wpp,cdf(GPR1_pd_t(19),wpp),'color',color_GPR,'linewidth',1.5)
plot(wpp,cdf('beta',wpp,a_beta(19),b_beta(19)),'color',color_MLE,'linewidth',1.5)
plot(wpp,cdf_aggr(19,:),'color',color_aggr,'linewidth',1.5,'linestyle','--')
plot(wind_power,p_settlement(:,19),'O','color',color_price)
legend('Agent1-ARMA','Agent2-GPR','Agent3-MLE','Aggregated in prediction market','Settled price in prediction market')
xlabel('Wind power (p.u.)')
ylabel('Cumulative Probability','rotation',0)
legend box off
box off
hold off
figure()
title ('hour 19')
%******plot price convergence in prediction market*****%
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot([0:17],squeeze(price2(5,19,:)))
hold on
plot([0:17],squeeze(price2(6,19,:)))
plot([0:17],squeeze(price2(7,19,:)))
xlabel('Trading Round')
ylabel(' Price in prediction market','rotation',0)
legend('If P_w < 0.5 p.u. ?','If P_w < 0.6 p.u. ?','If P_w < 0.7 p.u. ?')
legend box off
box off
%%
    real_value=[0.8385
    0.7067
    0.7281
    0.7509
    0.7215
    0.6910
    0.7900
    0.7588
    0.7211
    0.7588
    0.6763
    0.5855
    0.3746
    0.4040
    0.5941
    0.5726
    0.6496
    0.7182
    0.8194
    0.8062
    0.9525
    0.9600
    0.9527
    0.9458];
final_mat_bid=[ratio,real_value,bid_aggr',bid_agent1,bid_agent2,bid_agent3];
diff_aggr=[final_mat_bid(:,2)-final_mat_bid(:,3)];
diff_agent1=[final_mat_bid(:,2)-final_mat_bid(:,4)];
diff_agent2=[final_mat_bid(:,2)-final_mat_bid(:,5)];
diff_agent3=[final_mat_bid(:,2)-final_mat_bid(:,6)];
total_diff=[diff_aggr,diff_agent1,diff_agent2,diff_agent3];
%%
%%**** plot deviations ****
color_Arma=[0.78 0.08 0.52];
color_GPR=[0.12 0.56 1];
color_MLE=[0.13 0.55 0.13];
color_aggr=[0.1,0.1,0.44];
zer_h=zeros(1,length(diff_aggr));
figure()
box off
plot(diff_aggr,'color',color_aggr,'linewidth',2,'linestyle','--')
hold on
plot(diff_agent1,'color',color_Arma,'linewidth',1.5)
plot(diff_agent2,'color',color_GPR,'linewidth',1.5)
plot(diff_agent3,'color',color_MLE,'linewidth',1.5)
plot(zer_h,'linewidth',1)
hold off
xlabel('Hourly Period')
ylabel('Deviations (p.u.)', 'rotation',0)
legend('Aggregated in prediction market','Agent1-ARMA','Agent2-GPR','Agent3-MLE')
box off
legend boxoff 
index_error=sqrt(sum((total_diff.^2)))
% lamda_pos= xlsread(filename1);
% lamda_neg= xlsread(filename2);
%%
%************ plot bids and ratio ********%
color_Arma=[0.78 0.08 0.52];
color_GPR=[0.12 0.56 1];
color_MLE=[0.13 0.55 0.13];
color_aggr=[0.1,0.1,0.44];
figure()
%subplot(2,1,1)
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot(bid_aggr','color',color_aggr,'linewidth',2,'linestyle','--')
hold on
plot(bid_agent1,'color',color_Arma,'linewidth',1.5)
plot(bid_agent2,'color',color_GPR,'linewidth',1.5)
plot(bid_agent3,'color',color_MLE,'linewidth',1.5)
%plot(y_test(13:36),'k--','linewidth',2)
xlabel('Hourly Period')
ylabel('Bids (p.u.)', 'rotation',0)
legend('Aggregated in prediction market','Agent1-MLE','Agent2-GPR','Agent3-MLE','real generation')
legend box off
box off
%subplot(2,1,2)
%%
color_ratio=[0.98 0.5 0.45];
figure()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
stairs(ratio,'color',color_ratio,'linewidth',1.5)
xlabel('Hourly Period')
ylabel('\gamma - Percentile (%)','rotation',0)
box off

%%
%******hourly imbalance prices******%
filename2 = 'pos_after.xlsx';
filename3 = 'neg_after.xlsx';
lamda_pos_15min_after= xlsread(filename2);
lamda_neg_15min_after= xlsread(filename3);

lamda_pos_after2=zeros(24,1);
lamda_neg_after2=zeros(24,1);

for hp_after=1:24
    lamda_pos_after2(hp_after,1)=mean(lamda_pos_15min_after((hp_after-1)*4+1:hp_after*4));
end
for hn_after=1:24
    lamda_neg_after2(hn_after,1)=mean(lamda_neg_15min_after((hn_after-1)*4+1:hn_after*4));
end

typ2=lamda_pos_after2/max(lamda_pos_after2);
tyn2=lamda_neg_after2/max(lamda_neg_after2);
lamda_pos_after22=(1-typ2).*da_hourly_after;
lamda_neg_after22=(1+tyn2).*da_hourly_after;
lamda_pos_after2=da_hourly_after-lamda_pos_after22;
lamda_neg_after2=lamda_neg_after22-da_hourly_after;

for hi0=1:24
    if total_diff(hi0,1)>0
        imb1(hi0,1)=lamda_pos_after2(hi0,1).*total_diff(hi0,1);
    else
        imb1(hi0,1)=-lamda_neg_after2(hi0,1).*total_diff(hi0,1);
    end
end
imb_aggr=sum(imb1);
for hi1=1:24
    if total_diff(hi1,2)>0
        imb2(hi1,1)=lamda_pos_after2(hi1,1).*total_diff(hi1,2);
    else
        imb2(hi1,1)=-lamda_neg_after2(hi1,1).*total_diff(hi1,2);
    end
end
imb_arma=sum(imb2);

for hi2=1:24
    if total_diff(hi2,3)>0
        imb3(hi2,1)=lamda_pos_after2(hi2,1).*total_diff(hi2,3);                                                               
    else
        imb3(hi2,1)=-lamda_neg_after2(hi2,1).*total_diff(hi2,3);
    end
end
imb_gpr=sum(imb3);

for hi3=1:24
    if total_diff(hi3,4)>0
        imb4(hi3,1)=lamda_pos_after2(hi3,1).*total_diff(hi3,4);
    else
        imb4(hi3,1)=-lamda_neg_after2(hi3,1).*total_diff(hi3,4);
    end
end

imb_lme=sum(imb4);
imb_matrice1=[imb_aggr,imb_arma,imb_gpr,imb_lme]

%%
for hi0=1:24
    if total_diff(hi0,1)>0
        imb1(hi0,1)=lamda_pos(hi0,1).*total_diff(hi0,1);
    else
        imb1(hi0,1)=-lamda_neg(hi0,1).*total_diff(hi0,1);
    end
end
imb_aggr=sum(imb1);
for hi1=1:24
    if total_diff(hi1,2)>0
        imb2(hi1,1)=lamda_pos(hi1,1).*total_diff(hi1,2);
    else
        imb2(hi1,1)=-lamda_neg(hi1,1).*total_diff(hi1,2);
    end
end
imb_arma=sum(imb2);

for hi2=1:24
    if total_diff(hi2,3)>0
        imb3(hi2,1)=lamda_pos(hi2,1).*total_diff(hi2,3);                                                               
    else
        imb3(hi2,1)=-lamda_neg(hi2,1).*total_diff(hi2,3);
    end
end
imb_gpr=sum(imb3);

for hi3=1:24
    if total_diff(hi3,4)>0
        imb4(hi3,1)=lamda_pos(hi3,1).*total_diff(hi3,4);
    else
        imb4(hi3,1)=-lamda_neg(hi3,1).*total_diff(hi3,4);
    end
end

imb_lme=sum(imb4);
imb_matrice2=[imb_aggr,imb_arma,imb_gpr,imb_lme]

%%
%*************Fig 2: real generation**************%
Hourly_Wind_gen; 
%[y_pred_arma2,y_sd_arma2]=Arma2(); Plot Arma result **********
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
%subplot(2,2,1)
plot(y_train(end-2159:end),'k-','linewidth',2)
xticks([1:100:2160])
yticks([0:0.2:1])
xlim([1 2160])
xlabel('Training Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
box off
%%***** Fig. 2**** Three predictions*****
%%
color_Arma=[0.78 0.08 0.52];
color_GPR=[0.12 0.56 1];
color_MLE=[0.13 0.55 0.13];
color_stda=[1 0.89 1];
color_stdg=[0.69 0.88 0.9];
color_stdm=[0.6 0.98 0.6];
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 15; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
subplot(3,1,1)
ylim([0 1])
curve1=(y_pred_arma2-1*y_sd_arma2)';
curve2=(y_pred_arma2+1*y_sd_arma2)';
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1, fliplr(curve2)];
f=fill(x2', inBetween',color_stda,'LineStyle','none');
f.FaceAlpha=0.3;
hold on
plot(y_test(13:36),'k--','linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlim([1 24])
plot(y_pred_arma2,'color',color_Arma,'linewidth',2)
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
box off
legend('STD by Agent1-ARMA','Real generation','Mean by Agent1-ARMA')
legend box off
hold off
subplot(3,1,2)
%[y_pred2,y_sd2]=GPR2(); Plot GPR result **********
ylim([0 1])
curve1=(mu_GPR1-1*sigma_GPR1);
curve2=(mu_GPR1+1*sigma_GPR1);
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1, fliplr(curve2)];
f=fill(x2', inBetween',color_stdg,'LineStyle','none');
xticks([1:3:24])
yticks([0:0.2:1])
f.FaceAlpha=0.3;
hold on
plot(y_test(13:36),'k--','linewidth',2)
xticks([1:3:24])
yticks([0:0.2:1])
xlim([1 24])
plot(y_pred2,'color',color_GPR,'linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
legend('STD by Agent2-GPR','Real generation','Mean by Agent2-GPR')
legend box off
box off
subplot(3,1,3)
%[a_beta,b_beta]=mle22();****Plot MLE*****
mu_beta=a_beta./(a_beta+b_beta);
beta_sigma1=(a_beta.*b_beta)./((a_beta+b_beta+1).*(a_beta+b_beta).*(a_beta+b_beta));
beta_sigma=sqrt(beta_sigma1);
ylim([0 1])
curve1=(mu_beta-1*beta_sigma);
curve2=(mu_beta+1*beta_sigma);
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1, fliplr(curve2)];
f=fill(x2', inBetween',color_stdm,'LineStyle','none');
f.FaceAlpha=0.3;
yticks([0:0.2:1])
hold on
plot(y_test(13:36),'k--','linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlim([1 24])
plot(mu_beta,'color',color_MLE,'linewidth',2)
yticks([0:0.2:1])
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
legend('STD by Agent3-MLE','Real generation','Mean by Agent3-MLE')
legend box off
box off







