function [cnon,ceq] = agent_con(x,b,p,y1,z1)
[C,Q]=Cost_MM(b,p,x);
z2=z1+x; % z2 is the resulting asset of the agent after buying/selling x contract
y2=y1-C; % y2 is the resulting cash of the agent after buying/selling x contract
cnon(1)=-y2;
cnon(2)=-(z2+y2);
ceq=[];
end
