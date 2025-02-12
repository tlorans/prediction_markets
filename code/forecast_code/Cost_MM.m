function [C,Q]=Cost_MM(b,p,x)
%p is the current market price
%x is the number of contracts to purchase (x>0) or sell (x<0)
%b is the market maker parameter denoting the level of liquidity and
%sensitivity of the market price 
% C is the cost of buying/selling x contracts
% Q is the resulting market price after buying/selling x contracts 
C=b*log(p*(exp(x/b)-1)+1);
Q=1/(1+((1/p-1)/exp(x/b)));
end
