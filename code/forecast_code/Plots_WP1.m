function []=Plots_WP1(dayp)
%%
close all
dayp=163; %a number between 1 and 274
hp=10; %a number between 1 and 24
delta=12;
%*********** define colors ***************
color_da=[0.82 0.41 0.12];
color_pos=[1 0.41 0.71];
color_neg=[0.54 0.17 0.89];
color_Arma=[0.78 0.08 0.52];
color_GPR=[0.12 0.56 1];
color_MLE=[0.13 0.55 0.13];
color_stda=[1 0.89 1];
color_stdg=[0.69 0.88 0.9];
color_stdm=[0.6 0.98 0.6];
color_aggr=[0.1,0.1,0.44];
color_price=[1 0 0];
color_ratio=[0.98 0.5 0.45];
num_points=100;
wpp=linspace(0,1,num_points);
wind_power=0.1:0.1:0.9;
zer_h=zeros(1,24);
%***plot Fig .3 : electricity prices*********
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
stairs(Daily_in{5,dayp},'color',color_da,'linewidth',1.5)
hold on
stairs(Daily_in{6,dayp}{1},'color',color_pos,'linewidth',1.5)
stairs(Daily_in{6,dayp}{2},'color',color_neg,'linewidth',1.5)
ylabel('Electrcity Price (Euro/MWh)', 'rotation',0)
xlabel('Hour of the day')
legend ('day-ahead market clearing price','penalty price for positive deviation','penalty price for negative deviation')
box off
legend box off
%********hour=19****************%
figure()
title ('hour 19')
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot(wpp,cdf(Daily_out{1,dayp}.dis(hp+delta),wpp),'color',color_Arma,'linewidth',1.5)
hold on 
plot(wpp,cdf(Daily_out{2,dayp}.dis(hp+delta),wpp),'color',color_GPR,'linewidth',1.5)
plot(wpp,cdf('beta',wpp,Daily_out{3,dayp}.a(hp),Daily_out{3,dayp}.b(hp)),'color',color_MLE,'linewidth',1.5)
plot(wpp,Daily_out{4,dayp}.cdf(hp,:),'color',color_aggr,'linewidth',1.5,'linestyle','--')
plot(wind_power,Daily_out{4,dayp}.settled_price(:,hp),'O','color',color_price)
legend('Agent1-ARMA','Agent2-GPR','Agent3-MLE','Aggregated in prediction market','Settled price in prediction market')
xlabel('Wind power (p.u.)')
ylabel('Cumulative Probability','rotation',0)
legend box off
box off
hold off
figure()
title ('hour 10')
%******plot price convergence in prediction market*****%
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot([0:17],squeeze(Daily_out{4,dayp}.rounds_prices(5,hp,:)))
hold on
plot([0:17],squeeze(Daily_out{4,dayp}.rounds_prices(6,hp,:)))
plot([0:17],squeeze(Daily_out{4,dayp}.rounds_prices(7,hp,:)))
xlabel('Trading Round')
ylabel(' Price in prediction market','rotation',0)
legend('If P_w < 0.5 p.u. ?','If P_w < 0.6 p.u. ?','If P_w < 0.7 p.u. ?')
legend box off
box off

%%**** plot deviations ****
figure()
box off
plot(Daily_out{4,dayp}.deviation,'color',color_aggr,'linewidth',2,'linestyle','--')
hold on
plot(Daily_out{1,dayp}.deviation,'color',color_Arma,'linewidth',1.5)
plot(Daily_out{2,dayp}.deviation,'color',color_GPR,'linewidth',1.5)
plot(Daily_out{3,dayp}.deviation,'color',color_MLE,'linewidth',1.5)
plot(zer_h,'linewidth',1)
hold off
xlabel('Hourly Period')
ylabel('Deviations (p.u.)', 'rotation',0)
legend('Aggregated in prediction market','Agent1-ARMA','Agent2-GPR','Agent3-MLE')
box off
legend boxoff 

%************ plot bids and ratio ********%
figure()
%subplot(2,1,1)
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
plot(Daily_out{4,dayp}.bid,'color',color_aggr,'linewidth',2,'linestyle','--')
hold on
plot(Daily_out{1,dayp}.bid,'color',color_Arma,'linewidth',1.5)
plot(Daily_out{2,dayp}.bid,'color',color_GPR,'linewidth',1.5)
plot(Daily_out{3,dayp}.bid,'color',color_MLE,'linewidth',1.5)
%plot(y_test(13:36),'k--','linewidth',2)
xlabel('Hourly Period')
ylabel('Bids (p.u.)', 'rotation',0)
legend('Aggregated in prediction market','Agent1-MLE','Agent2-GPR','Agent3-MLE')
legend box off
box off
figure()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
stairs(Daily_in{7,dayp},'color',color_ratio,'linewidth',1.5)
xlabel('Hourly Period')
ylabel('\gamma - Percentile (%)','rotation',0)
box off
%*************Fig 2: real generation**************%
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
%subplot(2,2,1)
plot(Daily_in{2,dayp},'k-','linewidth',2)
xticks([1:100:2160])
yticks([0:0.2:1])
xlim([1 2160])
xlabel('Training Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
box off
%%***** Fig. 2**** Three predictions*****
figure ()
FigWidth = 8.89; % latex \linewdith of figure environment
FigHeight = 15; %Height/sqrt(2) which is a classicly nice proportion
set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
subplot(3,1,1)
ylim([0 1])
curve1=(Daily_out{1,dayp}.mu(13:36)-1*Daily_out{1,dayp}.sigm(13:36))';
curve2=(Daily_out{1,dayp}.mu(13:36)+1*Daily_out{1,dayp}.sigm(13:36))';
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1; fliplr(curve2)];
f=fill(x2', inBetween',color_stda,'LineStyle','none');
f.FaceAlpha=0.3;
hold on
plot(Daily_in{4,dayp}(13:36),'k--','linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlim([1 24])
plot(Daily_out{1,dayp}.mu(13:36),'color',color_Arma,'linewidth',2)
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
box off
legend('STD by Agent1-ARMA','Real generation','Mean by Agent1-ARMA')
legend box off
hold off
subplot(3,1,2)
%[y_pred2,y_sd2]=GPR2(); Plot GPR result **********
ylim([0 1])
curve1=(Daily_out{2,dayp}.mu(13:36)-1*Daily_out{2,dayp}.sigm(13:36));
curve2=(Daily_out{2,dayp}.mu(13:36)+1*Daily_out{2,dayp}.sigm(13:36));
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1; fliplr(curve2)];
f=fill(x2', inBetween',color_stdg,'LineStyle','none');
xticks([1:3:24])
yticks([0:0.2:1])
f.FaceAlpha=0.3;
hold on
plot(Daily_in{4,dayp}(13:36),'k--','linewidth',2)
xticks([1:3:24])
yticks([0:0.2:1])
xlim([1 24])
plot(Daily_out{2,dayp}.mu(13:36),'color',color_GPR,'linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
legend('STD by Agent2-GPR','Real generation','Mean by Agent2-GPR')
legend box off
box off
subplot(3,1,3)
%[a_beta,b_beta]=mle22();****Plot MLE*****
%mu_beta=a_beta./(a_beta+b_beta);
%beta_sigma1=(a_beta.*b_beta)./((a_beta+b_beta+1).*(a_beta+b_beta).*(a_beta+b_beta));
%beta_sigma=sqrt(beta_sigma1);
ylim([0 1])
curve1=(Daily_out{3,dayp}.mu-1*Daily_out{3,dayp}.sigm);
curve2=(Daily_out{3,dayp}.mu+1*Daily_out{3,dayp}.sigm);
x2 = [[1:24], fliplr([1:24])];
inBetween = [curve1; fliplr(curve2)];
f=fill(x2', inBetween',color_stdm,'LineStyle','none');
f.FaceAlpha=0.3;
yticks([0:0.2:1])
hold on
plot(Daily_in{4,dayp}(13:36),'k--','linewidth',2)
xticks([1:2:23])
yticks([0:0.2:1])
xlim([1 24])
plot(Daily_out{3,dayp}.mu,'color',color_MLE,'linewidth',2)
yticks([0:0.2:1])
xlabel('Hourly Period')
ylabel('Wind Power (p.u.)', 'rotation',0)
legend('STD by Agent3-MLE','Real generation','Mean by Agent3-MLE')
legend box off
box off
