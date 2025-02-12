function [WP,DA_price,IMB_price]=input_data()
%** Wind farms data one year-Australia 
filenameA = 'WP_5min_output_2018.xlsx'; 
xlRangeA = 'A2:A105121';
sheet1A = 'TARALGA1';
sheet3A = 'BALDHWF1';
sheet4A = 'WOOLNTH1';
subsetA1 = xlsread(filenameA,sheet1A,xlRangeA);
subsetA3 = xlsread(filenameA,sheet3A,xlRangeA);
load subsetA4_Saved.mat;
WP1_hourly=zeros(8760,1); 
WP3_hourly=zeros(8760,1);
WP4_hourly=zeros(8760,1);
for h=1:floor(length(subsetA1)/12)
    WP1_hourly(h,1)=mean(subsetA1((h-1)*12+1:h*12));
    WP3_hourly(h,1)=mean(subsetA3((h-1)*12+1:h*12));
end

for h=1:floor(length(subsetA4)/12)
    WP4_hourly(h,1)=mean(subsetA4((h-1)*12+1:h*12));
end


%%
%*****excluding negative values *********
WP1_hourly(WP1_hourly<0)=0;
WP3_hourly(WP3_hourly<0)=0;
WP4_hourly(WP4_hourly<0)=0;

field1='name';
value1={'Tarlaga','Bald Hills','Woolnorth'};
field2='capacity';
value2={107,107,140};
field3='location';
value3={'NSW','Victoria','Tasmania'};
field4='date';
value4={'2018'};
field5='hourly_output';
value5={WP1_hourly,WP3_hourly, WP4_hourly};
WP=struct(field1,value1,field3,value3,field4,value4);
%%
%** Electricity market data
%**NEM (National Electricity Market of Australia) day ahead prices
% filenameB = 'da_half_hour_price_2018.xlsx'; 
% xlRangeB = 'A2:L1489';
% sheet1B = 'NSW';
% sheet2B = 'VIC';
% sheet3B = 'TAS';
% subsetB1 = xlsread(filenameB,sheet1B,xlRangeB);
% subsetB2 = xlsread(filenameB,sheet2B,xlRangeB);
% subsetB3 = xlsread(filenameB,sheet3B,xlRangeB);

load Saved_subsetB.mat;

leng=12*length(subsetB1(:,1));
subsetB1=reshape(subsetB1,[leng,1]);
subsetB1(subsetB1==0)=[];
subsetB1(subsetB1<0)=0.01;

subsetB3=reshape(subsetB3,[leng,1]);
subsetB3(subsetB3==0)=[];
subsetB3(subsetB3<0)=0.01;


da_NSW=zeros(8760,1);
da_VIC=zeros(8760,1);
da_TAS=zeros(8760,1);

for h=1:8760
    da_NSW(h,1)=mean(subsetB1((h-1)*2+1:h*2));
    da_VIC(h,1)=mean(subsetB2((h-1)*2+1:h*2));
    da_TAS(h,1)=mean(subsetB3((h-1)*2+1:h*2));
    
end
field_d1='Day_ahead_Price';
value_d1={da_NSW,da_VIC,da_TAS};
DA_price=struct(field_d1,value_d1);


%%
%** Imbalance (deviation) prices 
filenameC = 'Imb_15min_price_2018.xls'; 
sheetsC=1:12;
xlRangeC_I='G3:G2978';
xlRangeC_D='H3:H2978';
MIP=xlsread(filenameC,sheetsC(1),xlRangeC_I);
MDP=xlsread(filenameC,sheetsC(1),xlRangeC_D);
for i=2:12
    file_temp_I=xlsread(filenameC,sheetsC(i),xlRangeC_I);
    MIP=[MIP;file_temp_I];
    file_temp_D=xlsread(filenameC,sheetsC(i),xlRangeC_D);
    MDP=[MDP;file_temp_D];
end
MDP(MDP<0)=0;
for h=1:length(MIP)/4
    hourly_MIP(h,1)=mean(MIP((h-1)*4+1:h*4));
    hourly_MDP(h,1)=mean(MDP((h-1)*4+1:h*4));
end

field1_imb='neg';
value1_imb={hourly_MDP};
field2_imb='pos';
value2_imb={hourly_MIP};
IMB_price=struct(field1_imb,value1_imb,field2_imb,value2_imb);








