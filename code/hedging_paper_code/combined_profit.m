function [profitw]=combined_profit(PM_low,PM_high,PM_Price,PM_NShares,Elec_bid,Imb_less,Imb_extra,Elec_daPrice,Cut_in,Vw_rated,Cut_off,Vwd,P_rated)
[gen]=pw_curve(Cut_in,Vw_rated,Cut_off,Vwd,P_rated);
Actual_revenue=zeros(size(Vwd,1),1);
[PM_Payoff]=Payoff_PMw(PM_low,PM_high,Vwd);

for i=1:size(gen,1)
    if gen(i)<= Elec_bid
        Actual_revenue(i)=Elec_bid*Elec_daPrice+(gen(i)-Elec_bid)*Imb_less;
    elseif gen(i)> Elec_bid
        Actual_revenue(i)=Elec_bid*Elec_daPrice+(gen(i)-Elec_bid)*Imb_extra;
    end
end

Actual_revenue=Elec_daPrice*gen;
profitw=Actual_revenue+(PM_NShares*(PM_Payoff-PM_Price));


