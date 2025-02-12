%%
close all
clear all
clc
%Input data
% pd = makedist('Weibul');
% pd.A=9;
% pd.B=1;
%Normal
%  pd = makedist('Normal');
%  pd.mu=6;
%  pd.sigma=0.3;
%Thid weibull works !
pd = makedist('weibull');
pd.a=8;
pd.b=10;
%Gamma
% pd = makedist('gamma');
% pd.a=6;
% pd.b=1;
 Vwd=random(pd,1e6,1);

%%
Cut_in=4;
Vw_rated=10;
Cut_off=30;
P_rated=1500;
%PM_low=6; % for normal
PM_low=median(pd); % for weibul
%PM_low=8.5;
%PM_low=6 % for gamma
%PM_high=pd.mu+2; %for normal
PM_high=pd.a+2;% for weibul and gamma
%PM_Price=0.5;
%PM_NShares=40;
Elec_bid=4000;
Imb_less=0.09;
Elec_daPrice=0.06;
Imb_extra=0.06;
%Wind Speed Plots and Wind Turbine Curves
[pw]=pw_curve(Cut_in,Vw_rated,Cut_off,Vwd,P_rated);
[x_range,Vwd_pdf,Vwd_cdf]=make_cdf(Vwd);

npricesamples=50;
%[nashpoint,Price_range,Win_EV,Loss_EV,exp_payoff]=find_nash(PM_low,PM_high,Vwd,npricesamples);
%nashpoint=0.5;
%Vwd_nash=(PM_high-PM_low)*nashpoint+PM_low;
[gen]=pw_curve(Cut_in,Vw_rated,Cut_off,Vwd,P_rated);
profitg=Elec_daPrice*gen;
[xg1_range,profitg_cm_pdf,profitg_cm_cdf]=make_cdf(profitg);
%Combined risks and profits vs price
%% calculating indifference and optimal pricing
PM_NShares1=10;
N=1*PM_NShares1;
b=0.1;
[min_v,index1_v]=min(abs(x_range-PM_low));
[min_v,index2_v]=min(abs(Vwd-PM_low));
[min_v,index1_g]=min(abs(xg1_range-Elec_daPrice*gen(index2_v)));
%p1g1=profitg_cm_cdf(index1_v)
Exp_R=trapz(xg1_range.*profitg_cm_pdf);
p1g2=profitg_cm_cdf(index1_g);
c=trapz([1:index1_g],xg1_range(1:index1_g).*profitg_cm_pdf(1:index1_g));
g=trapz([index1_g:length(xg1_range)],xg1_range(index1_g:length(xg1_range)).*profitg_cm_pdf(index1_g:length(xg1_range)));
%define p1
p11=Vwd_cdf(index1_v);
p1=cdf(pd,PM_low);
%indifference
%%
P_r_indiff=(log((p1*exp(-b*c)+(1-p1)*exp(-b*g))./(p1*exp(-b*N-b*c)+(1-p1)*exp(-b*g))))./(b*N);
%utility maximising
P_r_umax=(((1-p1)/p1)*exp(b*N-b*(g-c))+1).^(-1);
P_r_m_1=(((1-p1)/p1)*exp(b*N)+1).^(-1);
P_r_e1=(-log(p1*exp(-b*N)+(1-p1)))./(b*N);
%%
N=linspace(0.1,100,100);
P_r_indiff2=(log((p1*exp(-b*c)+(1-p1)*exp(-b*g))./(p1*exp(-b*N-b*c)+(1-p1)*exp(-b*g))))./(b*N);
%utility maximising
P_r_umax2=(((1-p1)/p1)*exp(b*N-b*(g-c))+1).^(-1);
P_r_m_12=(((1-p1)/p1)*exp(b*N)+1).^(-1);
P_r_e12=(-log(p1*exp(-b*N)+(1-p1)))./(b*N);
%%
set(0,'DefaultAxesFontSize',8); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','tex') %Allows us to use LaTeX maths notation
set(0, 'DefaultAxesFontName', 'times');
figure(2)
set(gcf, 'Units','centimeters')
set(gcf, 'Position',[0 0 8.89 3])
plot(N,P_r_indiff2,'-','linewidth',1.5,'color',[1,0.2,0.73]);
hold on
plot(N,P_r_umax2,'--','linewidth',1.5,'color',[1,0.2,0.73]);
plot(N,P_r_e12,'-','linewidth',1.5,'color',[0.3,0.53,1]);
plot(N,P_r_m_12,'--','linewidth',1.5,'color',[0.3,0.53,1]);
axis([0 100 0 0.8])
set(gca,'YTick',[0,0.2,0.5,0.8]) %Now impose sensible tickmark locations
set(gca,'YTickLAbel',{'0','0.2 $','0.5 $','0.8 $'})
set(gca, 'Units','centimeters')
set(gca, 'Position',[1.7 0.9 6.9 2])
ylabel({'Price'; 'of';'Shares'},'FontSiz',7.5) %Note cell matrices for line breaks
set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
xlabel('No. Shares','FontSiz',7.5) %Note use of unicode arrow for clarity
legend('Indifference pricing - WPP','Utility maximising - WPP','Indifference pricing - independent participant','Utility maximising - independent participant','position',[0.54,0.83,0.28,0.03])
hold off
box off
legend box off
%filename1 = ['fig1p', datestr(now,'mmmm dd yyyy hhMM'), '.png'] %Descriptive name timestamp and .png file format
%export_fig (filename1, '-m3', '-nocrop')

%%
PM_NShares1=10;
%Price_range=[P_r_m_1,p1,P_r_umax];
Price_range=[0.35,P_r_umax,0.75];
%profit
for i=1:size(Price_range,2)
    PM_Price=Price_range(i);
   [profitw]=combined_profit(PM_low,PM_high,PM_Price,PM_NShares1,Elec_bid,Imb_less,Imb_extra,Elec_daPrice,Cut_in,Vw_rated,Cut_off,Vwd,P_rated);
    profit_cm(:,i)=profitw;
end
for i=1:size(Price_range,2)
    [xp_range(:,i),profit_cm_pdf(:,i),profit_cm_cdf(:,i)]=make_cdf(profit_cm(:,i));
end
%Combined risks and profits vs number of Shares
PM_Price1=P_r_umax;
Share_range=linspace(1,10,3);

%profit
for i=1:size(Share_range,2)
    PM_NShares=Share_range(i);
   [profitw1]=combined_profit(PM_low,PM_high,PM_Price1,PM_NShares,Elec_bid,Imb_less,Imb_extra,Elec_daPrice,Cut_in,Vw_rated,Cut_off,Vwd,P_rated);
    profit1_cm(:,i)=profitw1;
end
for i=1:size(Share_range,2)
    [xp1_range(:,i),profit1_cm_pdf(:,i),profit1_cm_cdf(:,i)]=make_cdf(profit1_cm(:,i));
end


%%
%plot cdf of risk vs price
set(0,'DefaultAxesFontSize',8); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','tex') %Allows us to use LaTeX maths notation
set(0, 'DefaultAxesFontName', 'times');
figure (3)
set(gcf, 'Units','centimeters')
set(gcf,'Position',[0 0 8.89 3.3])%
plot(xg1_range,profitg_cm_cdf,'--','linewidth',1.5,'color',[0.5,0.5,0.5])
hold on
plot(xp_range(:,2),profit_cm_cdf(:,2),'-','linewidth',1.5,'color',[0,0.6,0.6])
plot(Elec_daPrice*gen(index2_v)-0.2,0.5,'*','color',[255,99,71]/255)
annotation('textbox','position',[0.174210631487279,0.36933334084352,0.249999997072986,0.279999992847443],'String','Increased revenues','FitBoxToText','on','FontName','times','FontSiz',6.5,'EdgeColor','none')
annotation('textbox','position',[0.714285714285714,0.543000004768372,0.139880949099149,0.239999995231629],'String','Reduced revenues','FitBoxToText','on','FontName','times','FontSiz',6.5,'EdgeColor','none')
annotation('doublearrow','position',[0.291666666666665,0.288,0.238095238095238,0.312])
annotation('doublearrow','position',[0.571428571428572,0.616,0.336309523809524,0.28])
%annotation('textbox',[0.3 0.2 .3 .4],'String',sprintf('Price of shares = %0.2f $',P_r_umax),'FitBoxToText','on','FontName','times','FontSiz',7,'EdgeColor','none')
axis([0 90 0 1])
set(gca,'YTick',[0,0.2,0.5,0.8,1]) %Now impose sensible tickmark locations
set(gca,'YTickLAbel',{'0','20 %','50 %','80 %','100%'})
set(gca,'XTick',[0:10:90]) %Now impose sensible tickmark locations
set(gca,'XTickLAbel',{'0','10 $','20 $','30 $','40 $','50 $','60 $', '70 $', '80 $', '90 $'})
set(gca, 'Units','centimeters')
set(gca, 'Position',[2.3 0.9 6.3 2.2])
ylabel({'Cumulative';'Probability'},'FontSiz',7.5) %Note cell matrices for line breaks
set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
xlabel('Revenue','FontSiz',7.5) %Note use of unicode arrow for clarity
legend('Unhedged','Hedged','Median','position',[0.299210631487279,0.68133334084352,0.249999997072986,0.279999992847443])
box off
legend box off
%%
set(0,'DefaultAxesFontSize',8); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','tex') %Allows us to use LaTeX maths notation
set(0, 'DefaultAxesFontName', 'times');
figure (4)
set(gcf, 'Units','centimeters')
set(gcf, 'Position',[0 0 8.89 3.3])%
plot(xp_range(:,1),profit_cm_cdf(:,1),'-','linewidth',1.2,'color',[0.8,0.8,1])
hold on
plot(xp_range(:,2),profit_cm_cdf(:,2),'-','linewidth',1.2,'color',[0.65,0.65,1])
plot(xp_range(:,3),profit_cm_cdf(:,3),'-','linewidth',1.2,'color',[0.35,0.4,1])
annotation('textbox',[0.3 0.3 .4 .5],'String',sprintf('No. Shares = %d',PM_NShares1),'FitBoxToText','on','FontName','times','FontSiz',7,'EdgeColor','none')
hold off
axis([0 90 0 1])
set(gca,'YTick',[0,0.2,0.5,0.8,1]) %Now impose sensible tickmark locations
set(gca,'YTickLAbel',{'0','20 %','50 %','80 %','100%'})
set(gca,'XTick',[0:10:90]) %Now impose sensible tickmark locations
set(gca,'XTickLAbel',{'0','10 $','20 $','30 $','40 $', '50 $', '60 $', '70 $', '80 $', '90 $'})
set(gca, 'Units','centimeters')
set(gca, 'Position',[2.3 0.9 6.3 2.2])
ylabel({'Cumulative';'Probability'},'FontSiz',7.5) %Note cell matrices for line breaks
set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
xlabel('Revenue','FontSiz',7.5) %Note use of unicode arrow for clarity
legend(sprintf('Price of shares = %0.2f $',Price_range(1)),sprintf('Price of shares = %0.2f $',Price_range(2)),sprintf('Price of shares = %0.2f $',Price_range(3)),'position',[0.77,0.53,0.01,0.03])
box off
legend box off
%filename5 = ['fig5p', datestr(now,'mmmm dd yyyy hhMM'), '.png'] %Descriptive name timestamp and .png file format
%export_fig (filename5, '-m3', '-nocrop')
set(0,'DefaultAxesFontSize',8); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','tex') %Allows us to use LaTeX maths notation
set(0, 'DefaultAxesFontName', 'times');
%%
figure (5) 
set(gcf, 'Units','centimeters')
set(gcf, 'Position',[0 0 8.89 3.3])%
plot(xp1_range(:,1),profit1_cm_cdf(:,1),'-','linewidth',1.2,'color',[1.00,0.8,0.8])
hold on
plot(xp1_range(:,2),profit1_cm_cdf(:,2),'-','linewidth',1.2,'color',[1.00,0.55,0.55])
plot(xp1_range(:,3),profit1_cm_cdf(:,3),'-','linewidth',1.2,'color',[1.00,0.3,0.3])
annotation('textbox',[0.28 0.3 .4 .6],'String',sprintf('Price of shares = %0.2f $',PM_Price1),'FitBoxToText','on','FontSiz',7,'FontName','times','EdgeColor','none')
hold off
axis([0 90 0 1])
set(gca,'YTick',[0,0.2,0.5,0.8,1]) %Now impose sensible tickmark locations
set(gca,'YTickLAbel',{'0','20 %','50 %','80 %','100%'})
set(gca,'XTick',[0:10:90]) %Now impose sensible tickmark locations
set(gca,'XTickLAbel',{'0','10 $','20 $','30 $','40 $', '50 $', '60 $', '70 $', '80 $', '90 $'})
set(gca, 'Units','centimeters')
set(gca, 'Position',[2.3 0.9 6.3 2.2])
ylabel({'Cumulative';'Probability'},'FontSiz',7.5) %Note cell matrices for line breaks
set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
xlabel('Revenue','FontSiz',7.5) %Note use of unicode arrow for clarity
legend(sprintf('No. Shares = %d',Share_range(1)),sprintf('No. Shares = %1.1f',Share_range(2)),sprintf('No. Shares = %d',Share_range(3)),'position',[0.75,0.53,0.01,0.03])
box off
legend box off
%filename4 = ['fig4p', datestr(now,'mmmm dd yyyy hhMM'), '.png'] %Descriptive name timestamp and .png file format
%export_fig (filename4, '-m3', '-nocrop')
%%
[U2]=p1*(-exp(-b*(N.*(1-P_r_umax2)+c)))+(1-p1)*(-exp(-b*(-N.*P_r_umax2+g)));
%[U1]=ones(1,length(U2))*(p1*(-exp(-b*c))+(1-p1)*(-exp(-b*g)));
[U3]=p1*(-exp(-b*(N.*(1-P_r_indiff2)+c)))+(1-p1)*(-exp(-b*(-N.*P_r_indiff2+g)));
set(0,'DefaultAxesFontSize',8); %Eight point Times is suitable typeface for an IEEE paper. Same as figure caption size
set(0,'DefaultFigureColor','w')
set(0,'defaulttextinterpreter','tex') %Allows us to use LaTeX maths notation
set(0, 'DefaultAxesFontName', 'times');
figure(6)
set(gcf, 'Units','centimeters')
set(gcf, 'Position',[0 0 8.89 2.5])% FigWidth = 8.89; % latex \linewdith of figure environment
plot(N,U2,'--','linewidth',1.5,'color',[1,0.2,0.73]);
hold on
plot(N,U3,'-','linewidth',1.5,'color',[1,0.2,0.73]);
axis([0 100 -0.2 -0.04])
set(gca, 'Units','centimeters')
set(gca, 'Position',[1.7 0.9 6.9 1.5])
ylabel({'Utility';'values'},'FontSiz',7.5) %Note cell matrices for line breaks
set(get(gca,'YLabel'),'Rotation',0, 'VerticalAlignment','middle', 'HorizontalAlignment','right') %Tidy it with right orientation (If all our vertical axes have the same internal offset all our axis labels will be neatly aligned
xlabel('No. Shares','FontSiz',7.5) %Note use of unicode arrow for clarity
legend('Utility maximising for WPP','Indifference pricing for WPP','position',[0.65,0.6,0.1,0.09])
hold off
box off
legend box off
%filename2 = ['fig2p', datestr(now,'mmmm dd yyyy hhMM'), '.png'] %Descriptive name timestamp and .png file format
%export_fig (filename2, '-m3', '-nocrop')
%%




