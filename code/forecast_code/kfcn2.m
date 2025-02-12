function y=kfcn2(XN,XM,theta)
% XN
% XM
% nparam=size(XN,2);
alfa=exp(theta);
% delta_v2=pdist2(XN(:,3),XM(:,3));
% delta_v1=pdist2(XN(:,2),XM(:,2));
delta_t=pdist2(XN(:,1),XM(:,1));
% kernel1_v1  =exp( -delta_v1.^2/2/alfa(1)^2 );
kernel1 =(1+delta_t.^2/2/alfa(1)/alfa(2)^2).^(-alfa(1));
kernel2  =exp( -delta_t.^2/2/alfa(3)^2 );
% kernel2_v2  =(1+delta_v2.^2/2/alfa(3)/alfa(4)^2).^(-alfa(3));
% kernel_t1   =(1+delta_t.^2/2/alfa(5)/alfa(6)^2).^(-alfa(5));
% kernel_t2   =exp( -delta_t .^2/2/alfa(7)^2 );
% kernel_t3   =exp(2*sin(delta_t/24*pi).^2/alfa(8)^2); %main time
% kernel_t4   =exp( -delta_t .^2/2/alfa(9)^2 );
% kernel_t5   =exp(2*sin(delta_t/24/7*pi).^2/alfa(10)^2); %main time

% kernel_e1=pdist2(XN(:,2),XM(:,2))<1e-15;
% kernel_e2=pdist2(XN(:,3),XM(:,3))<1e-15;

% size(kernel_e1);
% size(kernel1_v1);
% size(kernel_t2);
y=alfa(4)*kernel1+alfa(5)*kernel2 ;
% y=kernel1;
