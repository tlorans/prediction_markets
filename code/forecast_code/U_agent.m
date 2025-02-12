function [Exp_u]=U_agent(b,p,x,f,rho,z1,y1)
%p is the current market price
%x is the number of contracts to purchase (x>0) or sell (x<0)
%b is the market maker parameter denoting the level of liquidity and
%sensitivity of the market price 
% C is the cost of buying/selling x contracts
% Q is the resulting market price after buying/selling x contracts
%f is the belief of agent
[C,Q]=Cost_MM(b,p,x);
z2=z1+x; % z2 is the resulting asset of the agent after buying/selling x contract
y2=y1-C; % y2 is the resulting cash of the agent after buying/selling x contract
% agent utility function: U(W)=W^(1-rho)/(1-rho);
u1=((z2+y2).^(1-rho))/(1-rho);
u2=((y2).^(1-rho))/(1-rho);
% agent utility function: U(W)=-exp(-rhoW)
%u1=-exp((-rho)*(z2+y2));
%u2=-exp((-rho)*y2);
Exp_u=-(f*(u1)+(1-f)*(u2));
end
