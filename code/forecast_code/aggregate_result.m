function [cdf_aggr]=aggregate_result(p_settled,dfx)
% %****1)normal***%
% pd1 = makedist('Normal');
% pd1.mu=0.65;
% pd1.sigma=0.25;
% td1 = truncate(pd1,0,1);
% Exp_x1=mean(td1)
% x_td=linspace(0,1,10000);
% f_td1=@(x_td1) pdf(td1,x_td1);
% %***3)beta****%
% pd2 = makedist('beta');
% pd2.a=4;
% pd2.b=2;
% Exp_x2=mean(pd2)
% beta_mu=pd2.a/(pd2.a+pd2.b);
% beta_sigma=(pd2.a*pd2.b)/((pd2.a+pd2.b+1)*(pd2.a+pd2.b)*(pd2.a+pd2.b));
% f_pd2=@(x_td1) pdf(pd2,x_td1);
% %********5)rayleigh*****%
% pd3 = makedist('Rayleigh');
% pd3.B=0.6;
% td3 = truncate(pd3,0,1);
% Exp_x3=mean(td3)
% f_pd3=@(x_td1) pdf(td3,x_td1);
wind_power=0.1:0.1:0.9;

%***********
% p_settled = [0.0131    0.0290    0.0476    0.0751    0.1082    0.1483    0.1953    0.2493    0.3091    0.3765 0.4460    0.5177    0.5943    0.6698    0.7434    0.8171    0.8865];
% wind_power=0.1:0.05:0.9;


wind_power_adj=[0 wind_power 1];
p_settled_adj=[0 p_settled 1];
%weights=1./[ ...
    %0.0250    0.0229];
weights=1./dfx;   
weights_adj=[1000 weights 1000];

x=wind_power_adj;
y=p_settled_adj;
w=weights_adj;
func= @(b,x) betacdf(x,b(1),b(2));
% disp(x);
% disp(y);
% disp(w);

%modelfun = @(b,x)b(1) + b(2)*x(:,1).^b(3) + b(4)*x(:,2).^b(5);
%opts = statset('Display','iter','TolFun',1e-10);
opts = statset('Display','off','TolFun',1e-10);
try
    wnlm = fitnlm(x,y,func,[0.4 0.2],'Options',opts,'Weight',w);
catch
    wnlm = fitnlm(x,y,func,[0.3 0.7],'Options',opts,'Weight',w);
end

b1=table2array(wnlm.Coefficients(1,1));
b2=table2array(wnlm.Coefficients(2,1));
x_predict=linspace(0,1,100);
y_predict=betacdf(x_predict,b1,b2);

cdf_aggr=y_predict;
% y_agent1=cdf(td1,x_predict);
% y_agent2=cdf(pd2,x_predict);
% y_agent3=cdf(td3,x_predict);


% figure (5050)
% plot(x(2:end),y(2:end),'O')
% hold on
% plot(x_predict,y_predict,'r--')
% % plot(x_predict,y_agent1,'k-')
% % plot(x_predict,y_agent2,'g-')
% % plot(x_predict,y_agent3,'c-')
% hold off
% legend('PD market','WLSE')
% legend('PD market','WLSE','Agent1','Agent2','Agent3')
end
