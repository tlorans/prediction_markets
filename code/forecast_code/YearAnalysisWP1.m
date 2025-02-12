clear all
clc
close all
[WP,DA_price,IMB_price]=input_data();
which_WP=1;
hist_step=3*30*24; % historical data to train forecasting ~ 3 months
%days=(hist_step/24+1):(hist_step/24+5); % day that forecasting takes place for a look-ahead time equal to "forecast_step"
days=91:364;
%days=[91,96,126,210,250,301,364];
forecast_step=36;
useful_start=13;
useful_end=36;
counter=0;
num_meth=4;
bad_days=[];
Daily_cells=cell(length(days),1);
Daily_in=cell(7,length(days)); % rows: each input data, columns: each forecasting period
Daily_out=cell(num_meth,length(days)); % rows: forecsting methods, columns: each forecasting period
for i=1:length(days)
   counter=counter+1; 
   disp(days(counter))
   [C_input,C_output,bad_days]=Main_Execute(which_WP,WP,DA_price,IMB_price,days(i),forecast_step,useful_start,useful_end,hist_step,bad_days);
   Daily_cells(counter)={struct(("day"+counter),{C_input,C_output})};
   [Daily_in(:,counter),Daily_out(:,counter)]=Daily_cells{counter}.("day"+counter);
end


% hourly_count_min=37;
%hourly_count_max=8760-3*30*24-36;
% hourly_count_max=1+182*36;

% hourly_count_min=1+62*36;
% hourly_count_max=1+63*36;
% 
% num_days=(hourly_count_max-1)/36;



% count=1;
% for hourly_count=hourly_count_min:36:hourly_count_max
%     [C_input,C_output,bad_days]=Main_Execute(WP,DA_price,IMB_price,hourly_count,bad_days);
%     Daily_cells(count)={struct(("day"+count),{C_input,C_output})};
%     [Daily_in(:,count),Daily_out(:,count)]=Daily_cells{count}.("day"+count);
%     count=count+1
% end
% bad_days


% Daily_out{4,3}  % gives struct with the names of the fields
% Daily_out{4,3}.error % gives the values of the error for the 4th forecasting method in the third forecsting period
% Daily_in{1,2}  %give the value of the first set of iput data of the second forecasting period

%%
useful_step=useful_end-useful_start+1;
daily_deviations=zeros(useful_step,num_meth,length(days));
daily_error2=zeros(num_meth,length(days));
daily_error3=zeros(num_meth,length(days));
count=0;
for i=1:length(days)
    count=count+1;
    for meth=1:num_meth
 daily_error2(meth,count)=Daily_out{meth,count}.error(2); 
 daily_error3(meth,count)=Daily_out{meth,count}.error(3);  
 daily_deviations(:,meth,count)=Daily_out{meth,count}.deviation;  
    end
end


Daily_RMSE=100*sqrt(sum((daily_deviations.^2)))/sqrt(24);
total_deviation2=three_to_two(daily_deviations);
%total_RMSE2=100*sqrt(sum((total_deviation2.^2)))/sqrt(6576);
sum_dev=sum(total_deviation2);
sum_abs_dev=sum(abs(total_deviation2));
total_MBE=100*sum(total_deviation2)/(length(days)*24);
total_MAE=100*sum(abs(total_deviation2))/(length(days)*24);
total_MSE=100*sum((total_deviation2.^2))/(length(days)*24);
total_RMSE2=100*sqrt(sum((total_deviation2.^2))/(length(days)*24));
total_error2=sum(daily_error2');
total_error3=sum(daily_error3');
%%
sup_dev=0;
for j=1:length(days)
    for k=1:24
   if Daily_in{6,j}{2}(k)>Daily_in{6,j}{1}(k)
       sup_dev=sup_dev+1;
   end
    end
end
%%
count_daily_RMSE=0;
for j=1:length(days)
   if min([Daily_RMSE(1,:,j)])==Daily_RMSE(1,4,j)
       count_daily_RMSE=count_daily_RMSE+1;
    end
end


%%
save myWP1



