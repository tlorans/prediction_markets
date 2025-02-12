function [C_input,C_output,bad_days]=Main_Execute(which_WP,WP,DA_price,IMB_price,which_day,forecast_step,useful_start,useful_end,hist_step,bad_days)
%clc
%clear all
%close all
%%
%******Part 1: Load Data  
%first_h_day=3430-57-107; % Date analysed in the first paper 
first_h_day=which_day*24-(useful_start-1)-hist_step;
WP_cap=WP(which_WP).capacity;
%find(isnan(X)==1)
%%
hist_imb=which_day*24-4*24; %historical mean duration for forecasting : 4 days
da_hourly_after=DA_price(which_WP).Day_ahead_Price((which_day*24):((which_day*24)+23));
hourly_MDP=IMB_price.neg((hist_imb:(which_day*24)),1);
hourly_MID=IMB_price.pos((hist_imb:(which_day*24)),1);
for h=1:24
    daily_MID(h,1)=mean(hourly_MID(h:24:(end-24+h)));
    daily_MDP(h,1)=mean(hourly_MDP(h:24:(end-24+h)));
end
lamda_pos_after=daily_MID;
lamda_neg_after=daily_MDP;
typ=lamda_pos_after/max(lamda_pos_after);
tyn=lamda_neg_after/max(lamda_neg_after);
lamda_pos_after=(1-typ).*da_hourly_after;
lamda_neg_after=(1+tyn).*da_hourly_after;
lamda_pos=da_hourly_after-lamda_pos_after;
lamda_neg=lamda_neg_after-da_hourly_after;
ratio=lamda_pos./(lamda_pos+lamda_neg);
%%
%***************imbalance prices at the target forecasting day ****
lamda_pos_after2=IMB_price.neg((which_day*24):((which_day*24)+23));
lamda_neg_after2=IMB_price.pos((which_day*24):((which_day*24)+23));
typ2=lamda_pos_after2/max(lamda_pos_after2);
tyn2=lamda_neg_after2/max(lamda_neg_after2);
lamda_pos_after22=(1-typ2).*da_hourly_after;
lamda_neg_after22=(1+tyn2).*da_hourly_after;
lamda_pos_after2=da_hourly_after-lamda_pos_after22;
lamda_neg_after2=lamda_neg_after22-da_hourly_after;
%*********cell array containing historical and the operaring day imbalance prices
imbalance_cell={lamda_pos,lamda_neg,lamda_pos_after2,lamda_neg_after2};
%%
%******Part 2: Train dataset
[x_train,y_train,x_test,y_test]=train_test_set(which_WP,first_h_day,WP,forecast_step,hist_step);
%***** Part 3: Agents Predict WP Outputs
[y_pred_arma2,y_sd_arma2,bad_days_temp]=Arma2(y_train,which_day,forecast_step);
bad_days=[bad_days;bad_days_temp];
[y_pred2,y_sd2]=GPR2(x_train,y_train,x_test,hist_step);
%[a_beta,b_beta]=mle22(which_WP,first_h_day,WP,forecast_step,hist_step);
[a_beta,b_beta]=mle22(which_WP,first_h_day,WP,hist_step);
mu_beta=a_beta./(a_beta+b_beta);
beta_sigma1=(a_beta.*b_beta)./((a_beta+b_beta+1).*(a_beta+b_beta).*(a_beta+b_beta));
beta_sigma=sqrt(beta_sigma1);
%******** Part 4: bid values
useful_step=useful_end-useful_start+1;
bid_agent1=zeros(useful_step,1);
bid_agent2=zeros(useful_step,1);
bid_agent3=zeros(useful_step,1);
num_points=100;
wpp=linspace(0,1,num_points);
wind_power=0.1:0.1:0.9;
b_1=20;
p_settlement=zeros(length(wind_power),useful_step);
for hr=useful_start:useful_end
    which_val=ratio(hr-useful_start+1);    
    %cdf_Arma1 = cdf('Normal',wpp,y_pred_arma2(hr),y_sd_arma2(hr));
    Arma1_pd(hr) = makedist('Normal','mu',y_pred_arma2(hr),'sigma',y_sd_arma2(hr));
    Arma1_pd_t(hr)= truncate(Arma1_pd(hr),0,1);
    mu_Arma1(hr)=mean(Arma1_pd_t(hr));
    sigma_Arma1(hr)=std(Arma1_pd_t(hr));
    cdf_Arma1= cdf(Arma1_pd_t(hr),wpp);
    distance_arma=(cdf_Arma1-which_val).^2;
    [~,find_val_arma]=min(distance_arma);
    bid_agent1(hr-useful_start+1)=wpp(find_val_arma);
    GPR1_pd(hr) = makedist('Normal','mu',y_pred2(hr),'sigma',y_sd2(hr));
    GPR1_pd_t(hr)= truncate(GPR1_pd(hr),0,1);
    mu_GPR1(hr)=mean(GPR1_pd_t(hr));
    sigma_GPR1(hr)=std(GPR1_pd_t(hr));
    %cdf_GPR1= cdf('Normal',wpp,y_pred2(hr),y_sd2(hr));
    %cdf_GPR1= cdf('Normal',wpp,mu_GPR1(hr),sigma_GPR1(hr));
    cdf_GPR1= cdf(GPR1_pd_t(hr),wpp);
    distance_GPR=(cdf_GPR1-which_val).^2;
    [~,find_val_GPR]=min(distance_GPR);
    bid_agent2(hr-useful_start+1)=wpp(find_val_GPR);
    cdf_mle1= cdf('beta',wpp,a_beta(hr-useful_start+1),b_beta(hr-useful_start+1));
    distance_mle=(cdf_mle1-which_val).^2;
    [~,find_val_mle]=min(distance_mle);
    bid_agent3(hr-useful_start+1)=wpp(find_val_mle);
    % Aggregate cdf
    for kb=1:length(wind_power)
        which_vall=wind_power(kb);
        distance=(wpp-which_vall).^2;
        [~,find_vall]=min(distance);
        agent_1(kb)=cdf_Arma1(find_vall);
        agent_2(kb)=cdf_GPR1(find_vall);
        agent_3(kb)=cdf_mle1(find_vall);
        [p_settlement(kb,hr-useful_start+1),price2(kb,hr-useful_start+1,:)]=traders_models(agent_1(kb),agent_2(kb),agent_3(kb),b_1);
        p_settled(kb)=p_settlement(kb,hr-useful_start+1);
        dfx(kb)=(exp(0/b_1).*(1/p_settled(kb) - 1))./(b_1.*(exp(0./b_1).*(1/p_settled(kb) - 1) + 1).^2);
        weight_h(kb,hr-useful_start+1)=dfx(kb);
    end
    cdf_aggr(hr-useful_start+1,:)=aggregate_result(p_settled,dfx);
    distance_aggr=(cdf_aggr(hr-useful_start+1,:)-which_val).^2;
    [~,find_val_aggr]=min(distance_aggr);
    bid_aggr(hr-useful_start+1)=wpp(find_val_aggr);
end

%*********gather results of the forecasting methods for the evluation (assessment)
%%
method1='arma';
method2='gpr';
method3='mle';
method4='aggr'; %prediction market
bid_value1=bid_agent1;
bid_value2=bid_agent2;
bid_value3=bid_agent3;
bid_value4=bid_aggr;
methods_bids=struct(method1,bid_value1,method2,bid_value2,method3,bid_value3,method4,bid_value4);
[index_error,total_diff]=assessment(WP_cap,y_test,methods_bids,imbalance_cell,useful_start,useful_end);
%********* All results for one day (24 hours) in each forecasting method ********
arma_s=struct('mu',mu_Arma1,'sigm',sigma_Arma1,'dis',Arma1_pd_t,'bid',bid_agent1,'deviation',total_diff(:,2),'error',index_error.arma);
gpr_s=struct('mu',mu_GPR1','sigm',sigma_GPR1','dis',GPR1_pd_t,'bid',bid_agent2,'deviation',total_diff(:,3),'error',index_error.gpr);
mle_s=struct('a',a_beta','b',b_beta','mu',mu_beta','sigm',beta_sigma','bid',bid_agent3,'deviation',total_diff(:,4),'error',index_error.mle);
aggr_s=struct('settled_price',p_settlement,'rounds_prices',price2,'weights',weight_h,'cdf',cdf_aggr,'bid',bid_aggr','deviation',total_diff(:,1),'error',index_error.aggr);
%********* cell arrays of the final input output *************%
C_input={x_train;y_train;x_test;y_test;da_hourly_after;imbalance_cell;ratio};
C_output={arma_s;gpr_s;mle_s;aggr_s};











