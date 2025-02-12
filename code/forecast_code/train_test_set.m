
function [x_train,y_train,x_test,y_test]=train_test_set(which_WP,first_h_day,WP,forecast_step,hist_step)
%a1=first_h_day-forecast_step;
a1=first_h_day;
%a2=first_h_day+hist_step-forecast_step; % three months training period
a2=first_h_day+hist_step; % three months training period
b1=a2;
b2=b1+forecast_step;
WP_cap=WP(which_WP).capacity;
WP_hourly=WP(which_WP).hourly_output;
x_train=(1:1:(a2-a1))';
y_train=WP_hourly(a1:a2-1)/WP_cap;
x_test=((a2-a1)+1:1:(a2-a1)+forecast_step)';
y_test=WP_hourly(b1:b2-1)/WP_cap;



