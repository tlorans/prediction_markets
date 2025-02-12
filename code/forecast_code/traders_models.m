function [p_settlement,price2]=traders_models(agent_1,agent_2,agent_3,b_1)
close all
%clear all
clc
%*****sample subjective beliefs of the agents***%%%
%****1)normal***%
% pd1 = makedist('Normal');
% pd1.mu=0.65;
% pd1.sigma=0.25;
% td1 = truncate(pd1,0,1);
% Exp_x1=mean(td1)
% x_td=linspace(0,1,10000);
% f_td1=@(x_td1) pdf(td1,x_td1);
%***3)beta****%
% pd2 = makedist('beta');
% pd2.a=4;
% pd2.b=2;
% Exp_x2=mean(pd2)
% beta_mu=pd2.a/(pd2.a+pd2.b);
% beta_sigma=(pd2.a*pd2.b)/((pd2.a+pd2.b+1)*(pd2.a+pd2.b)*(pd2.a+pd2.b));
% f_pd2=@(x_td1) pdf(pd2,x_td1);
%******3)log normal*******%
%pd3 = makedist('Lognormal');
%pd3.mu=log((0.58.^2)/sqrt(0.58^2+0.06^2));
%pd3.sigma=sqrt(log(1+0.06^2/0.58.^2));
%Exp_x3=mean(pd3)
%f_pd3=@(x_td1) pdf(pd3,x_td1);
%*********4)*************************%
% %pd4 = makedist('Normal');
% %pd4.mu=0.58;
% pd4.sigma=0.06;
% td4 = truncate(pd4,0,1);
% Exp_x4=mean(td4)
% f_pd4=@(x_td4) pdf(td4,x_td4);
%********5)rayleigh*****%
% pd3 = makedist('Rayleigh');
% pd3.B=0.6;
% td3 = truncate(pd3,0,1);
% Exp_x3=mean(td3)
% f_pd3=@(x_td1) pdf(td3,x_td1);
%***********************************%




%bid=bid_1;
%agent1=integral(f_td1,0,bid)
%agent2=integral(f_pd2,0,bid)
%agent3=integral(f_pd3,0,bid)
%agent4=integral(f_pd4,0,bid)
%agent5=integral(f_pd5,0,bid)
% agent1=cdf(td1,bid);
% agent2=cdf(pd2,bid);
% agent3=cdf(td3,bid);
agent1=agent_1;
agent2=agent_2;
agent3=agent_3;
%************************************
n=3; % number of agents (participants/traders)
pi=[agent1,agent2,agent3]; % agents beliefs of the outcome (subjective probability)
%pi=rand(1,n);%agents beliefs uniformly distributed  in [0,1]
yi0=[50,50,50];%agents initial cash (budget)
%yi0=10+10.*rand(1,n); %initial budgets (cash) uniformly distributed  in [10,20]
zi0=zeros(1,n);%agents assets
rho=2;
yit=yi0;
zit=zi0;%assets of each agent
rit=zi0;%number of contracts traded by each agent
%p0=min(min(agent1,agent2),agent3); % market maker initial belief
p0=0.5;
p00=p0;
t=17; % number of trading periods
b=b_1; % LMSR liqudity Parameter
qt=0;
A=[]; 
b1=[];
Aeq=[];
beq=[];
x0=qt;
lb=[];
ub=[];
[d_series]=trading_order(t,n); 
price=zeros(t,1);%market price in each tradind period 
quantity=zeros(t,1);%number of contracts in market in each tradind period 
qx=zeros(t,n);%number of contracts traded by each agent in each period
fi=ones(t,1)*pi;%beliefs of agents in each period
flag2=zeros(t,n);
alfa=0;%learning rate for revising agents' beliefs according to market price
kt=zeros(t,1)*[1:n]; %for denoting in which rounds, each agent accesses the market
m_last=zeros(2,n); % for recording last time traded by each agent
count=0;
for k=1:t
   j=d_series(k); %specify which agent access the market at this trading period (round)
   %j=randi([1 n],1,1);
   if k==1 %record the rounds of trading of each agent
   kt(k,j)=k;
   end
   if k>1
   kt(k,:)=kt(k-1,:);
   kt(k,j)=k;
   end
 m_last(1,j)=m_last(2,j);  
 m_last(2,j)=k; 
 %revise belief
 if k>1
 p2=price(k-1,1); % current market price
 else
     p2=p00;
 end
  if m_last(1,j)~=0
  p1=price(m_last(1,j),1);
  f1=fi(m_last(1,j),j);
  else
     p1=p00;
     f1=pi(j); 
  end
 if k>1
 fi(k,:)=fi(k-1,:);
 end
 if k==1
    fi(k,:)= pi;
 end
 [f2,flg]=revise_belief(p1,p2,f1,alfa);
 fi(k,j)=f2;
 flag2(k,j)=flg; %if the belief of trader j at period k is revised
 f=fi(k,j);
 %f=pi(j);
 p=p0;
 z1=zit(j);
 y1=yit(j);
       %[x,fval] = fmincon(@(x)Exp_agent(b,p,x,f),x0,A,b1,Aeq,beq,lb,ub,@(x)agent_con(x,b,p,y1,z1));
       opt = optimset('Display','off');
      [x,fval] = fmincon(@(x)U_agent(b,p,x,f,rho,z1,y1),x0,A,b1,Aeq,beq,lb,ub,@(x)agent_con(x,b,p,y1,z1),opt);
       rit(j)=x;
       qx(k,j)=rit(j);           
       zit(j)=zit(j)+rit(j);
       [C,Q]=Cost_MM(b,p,rit(j));
       yit(j)=yit(j)-C;
       qt=qt+rit(j);
       p0=Q;
       price(k,1)=Q;
       quantity(k,1)=qt;         
end
%price(t,1)
p_settlement=price(t,1);
y00=b*log(1/min(p00,1-p00));
yit2=[y00,yi0];
y_bar=sum(yit2);
price2=[p00;price];
p002=[p00,pi];
p_bar=(1/y_bar)*sum(yit2.*p002);%market price converge to budget weighted average of traders beliefs

% figure()
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot([0:t],price2)
% %legend(sprintf('N=%3.2f',0'),sprintf('N^*_m=%3.2f',N(1)),sprintf('N^*=%3.2f',N(2)),'location','best')
% xlabel('Trading Round')
% %ylabel('Price','rotation',0)
% ylabel('Price')
% box off
%legend boxoff 


% xagg=[0
% 0.15
% 0.25
% 0.35
% 0.45
% 0.55
% 0.65
% 0.75
% 0.85
% 0.95
% 1];
% yagg=[0
% 0.1171
% 0.1907
% 0.274
% 0.3634
% 0.457
% 0.5536
% 0.6535
% 0.7588
% 0.8766
% 1
% ];

%vq=interp1(xagg,yagg,x_td,'spline');

% figure(2)
% FigWidth = 8.89; % latex \linewdith of figure environment
% FigHeight = 5; %Height/sqrt(2) which is a classicly nice proportion
% set(gcf, 'Units','centimeters', 'Position',[0 0 FigWidth FigHeight]) %Impose print size on figure
% set(gcf,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[FigWidth FigHeight])
% plot(x_td, cdf(td1,x_td));
% hold on
% plot(x_td, cdf(pd2,x_td));
% plot(x_td, cdf(td3,x_td));
% %plot(x_td, pdf(td4,x_td));
% %plot(x_td, pdf(pd5,x_td));
% %plot(vq,x_td);
% legend('agent 1','agent 2', 'agent3')
% hold off
% xlabel('Wind Power (p.u.)')
% ylabel('Cumulative Probability')
% box off
% legend box off 

% agent1
% agent2
% agent3
end



