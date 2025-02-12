function [index_error,total_diff]=assessment(WP_cap,y_test,methods_bids,imbalance_cell,useful_start,useful_end)
%%
%real_value=y_test(13:36);
real_value=y_test(useful_start:useful_end);
%imbalance_cell={lamda_pos,lamda_neg,lamda_pos_after2,lamda_neg_after2};
lamda_pos=imbalance_cell{1};
lamda_neg=imbalance_cell{2};
lamda_pos_after2=imbalance_cell{3};
lamda_neg_after2=imbalance_cell{4};

method1='arma';
method2='gpr';
method3='mle';
method4='aggr'; %prediction market

bid_aggr=methods_bids.aggr;
bid_agent1=methods_bids.arma;
bid_agent2=methods_bids.gpr;
bid_agent3=methods_bids.mle;


final_mat_bid=[real_value,bid_aggr',bid_agent1,bid_agent2,bid_agent3];
diff_aggr=[final_mat_bid(:,1)-final_mat_bid(:,2)];
diff_agent1=[final_mat_bid(:,1)-final_mat_bid(:,3)];
diff_agent2=[final_mat_bid(:,1)-final_mat_bid(:,4)];
diff_agent3=[final_mat_bid(:,1)-final_mat_bid(:,5)];
total_diff=[diff_aggr,diff_agent1,diff_agent2,diff_agent3];
%*****index1: RMSE with respect to bids in each method *********
index_error=sqrt(sum((total_diff.^2)));
index1=100*index_error'/sqrt(24);
%*****index2: imbalance costs ********
%************assess based on actual prices
for hi0=1:24
    if total_diff(hi0,1)>0
        imb1(hi0,1)=lamda_pos_after2(hi0,1).*total_diff(hi0,1);
    else
        imb1(hi0,1)=-lamda_neg_after2(hi0,1).*total_diff(hi0,1);
    end
end
imb_aggr=sum(imb1);
for hi1=1:24
    if total_diff(hi1,2)>0
        imb2(hi1,1)=lamda_pos_after2(hi1,1).*total_diff(hi1,2);
    else
        imb2(hi1,1)=-lamda_neg_after2(hi1,1).*total_diff(hi1,2);
    end
end
imb_arma=sum(imb2);

for hi2=1:24
    if total_diff(hi2,3)>0
        imb3(hi2,1)=lamda_pos_after2(hi2,1).*total_diff(hi2,3);                                                               
    else
        imb3(hi2,1)=-lamda_neg_after2(hi2,1).*total_diff(hi2,3);
    end
end
imb_gpr=sum(imb3);

for hi3=1:24
    if total_diff(hi3,4)>0
        imb4(hi3,1)=lamda_pos_after2(hi3,1).*total_diff(hi3,4);
    else
        imb4(hi3,1)=-lamda_neg_after2(hi3,1).*total_diff(hi3,4);
    end
end

imb_lme=sum(imb4);
imb_matrice1=[imb_aggr,imb_arma,imb_gpr,imb_lme];
index2=WP_cap*imb_matrice1';
%******* assess based on estimated prices
for hi0=1:24
    if total_diff(hi0,1)>0
        imb1(hi0,1)=lamda_pos(hi0,1).*total_diff(hi0,1);
    else
        imb1(hi0,1)=-lamda_neg(hi0,1).*total_diff(hi0,1);
    end
end
imb_aggr=sum(imb1);
for hi1=1:24
    if total_diff(hi1,2)>0
        imb2(hi1,1)=lamda_pos(hi1,1).*total_diff(hi1,2);
    else
        imb2(hi1,1)=-lamda_neg(hi1,1).*total_diff(hi1,2);
    end
end
imb_arma=sum(imb2);

for hi2=1:24
    if total_diff(hi2,3)>0
        imb3(hi2,1)=lamda_pos(hi2,1).*total_diff(hi2,3);                                                               
    else
        imb3(hi2,1)=-lamda_neg(hi2,1).*total_diff(hi2,3);
    end
end
imb_gpr=sum(imb3);

for hi3=1:24
    if total_diff(hi3,4)>0
        imb4(hi3,1)=lamda_pos(hi3,1).*total_diff(hi3,4);
    else
        imb4(hi3,1)=-lamda_neg(hi3,1).*total_diff(hi3,4);
    end
end

imb_lme=sum(imb4);
imb_matrice2=[imb_aggr,imb_arma,imb_gpr,imb_lme];
index3=WP_cap*imb_matrice2';
method1='aggr';
method2='arma';
method3='gpr';
method4='mle';
v1=[index1(1);index2(1);index3(1)];
v2=[index1(2);index2(2);index3(2)];
v3=[index1(3);index2(3);index3(3)];
v4=[index1(4);index2(4);index3(4)];
index_error=struct(method1,v1,method2,v2,method3,v3,method4,v4);



