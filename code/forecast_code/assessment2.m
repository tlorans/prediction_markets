function [index_error,total_diff]=assessment2(WP_cap,y_test,bid_simple,imbalance_cell,useful_start,useful_end)
%%
real_value=y_test(useful_start:useful_end);
%imbalance_cell={lamda_pos,lamda_neg,lamda_pos_after2,lamda_neg_after2};
lamda_pos=imbalance_cell{1};
lamda_neg=imbalance_cell{2};
lamda_pos_after2=imbalance_cell{3};
lamda_neg_after2=imbalance_cell{4};


bid_ens=bid_simple;
final_mat_bid=[real_value,bid_ens];
diff_ens=[final_mat_bid(:,1)-final_mat_bid(:,2)];
total_diff=diff_ens;
%*****index1: RMSE with respect to bids in each method *********
index_error=sqrt(sum((total_diff.^2)));
index1=100*index_error'/sqrt(24);
%*****index2: imbalance costs ********
%************assess based on actual prices

for hi3=1:24
    if total_diff(hi3,1)>0
        imb4(hi3,1)=lamda_pos_after2(hi3,1).*total_diff(hi3,1);
    else
        imb4(hi3,1)=-lamda_neg_after2(hi3,1).*total_diff(hi3,1);
    end
end

imb_lme=sum(imb4);
imb_matrice1=imb_lme;
index2=WP_cap*imb_matrice1';
%******* assess based on estimated prices

for hi3=1:24
    if total_diff(hi3,1)>0
        imb4(hi3,1)=lamda_pos(hi3,1).*total_diff(hi3,1);
    else
        imb4(hi3,1)=-lamda_neg(hi3,1).*total_diff(hi3,1);
    end
end

imb_lme=sum(imb4);
imb_matrice2=imb_lme;
index3=WP_cap*imb_matrice2';
index_error=[index1;index2;index3];


