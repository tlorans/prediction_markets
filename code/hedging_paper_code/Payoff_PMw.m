function [PM_Payoff]=Payoff_PMw(PM_low,PM_high,Vwd)
PM_Payoff=zeros(size(Vwd,1),1);
for i=1:size(Vwd,1)
    if Vwd(i)< PM_low
        PM_Payoff(i)=1;
    else
       PM_Payoff(i)=0; 
    end
end

%PM_Payoff=1-PM_Payoff;

                                                                                              