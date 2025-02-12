function [pw]=pw_curve(Cut_in,Vw_rated,Cut_off,Vwd,P_rated)
pw=zeros(size(Vwd,1),1);
for i=1:size(Vwd,1)
    if Vwd(i)<Cut_in
        pw(i)=0;
    elseif (Vwd(i)>=Cut_in)&&(Vwd(i)<Vw_rated)
        pw(i)=P_rated/(Vw_rated^3-Cut_in^3)*(Vwd(i)^3-Cut_in^3);
    elseif (Vwd(i)>=Vw_rated)&&(Vwd(i)<Cut_off)
        pw(i)=P_rated;
    elseif Vwd(i)>Cut_off
        pw(i)=0;
    end
end
        
