function [y_pred2,y_sd2]=GPR2(x_train,y_train,x_test,hist_step)
%range1=1560:2160;
%%
range1=(hist_step-600):(hist_step);
%%
% theta0=log(10*[0.3 0.3 0.3 0.3 0.3 500 500]);
theta0=log(50*[0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 500 500 500]);
%theta0=log(10*[0.3 0.3 0.3 0.3 0.3 0.3 0.3 0.3 50 50 200]);
sigma0=0.2;
%sigma0=0.02;
fun=@kfcn7;
gprMdl = fitrgp(x_train(range1),y_train(range1),'KernelFunction',fun,'KernelParameters',theta0,'sigma',sigma0);
calc_LL=gprMdl.LogLikelihood;
% param1=gprMdl1.HyperparameterOptimizationResults;
kernelparam=gprMdl.KernelInformation.KernelParameters;
L = resubLoss(gprMdl);
[y_pred,y_sd,~] = predict(gprMdl,x_test);
[y_pred_train,y_sd_train,~] = predict(gprMdl,x_train(range1));
y_pred2=y_pred;
y_sd2=y_sd;
end



