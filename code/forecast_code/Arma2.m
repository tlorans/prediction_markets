function [y_pred_arma2,y_sd_arma2,bad_days]=Arma2(y_train,which_day,forecast_step)
num_y=72;
num_e=36;
try
    ARMA_mdl1 = arima('ARLags',1:num_y,'MALags',1:num_e,'Constant',0);
    EstMdl = estimate(ARMA_mdl1,y_train,'Display','off');
    [y_pred_arma,y_sd_arma] = forecast(EstMdl,forecast_step,y_train(end-5*num_y:end));
    bad_days=[];
catch
    y_pred_arma=rand(forecast_step,1);
    y_sd_arma=rand(forecast_step,1);
    bad_days=which_day;
end
y_pred_arma2=y_pred_arma;
y_sd_arma2=y_sd_arma;
end

