%two dimensional search
clear;%clear all data
clc;%clear the screen

Q_h_vector=10:0.1:60;
num1_to_search=length(Q_h_vector);
Q_l_vector=10:0.1:60;
num2_to_search=length(Q_l_vector);

res_mat_1=zeros(num1_to_search,num2_to_search);%diff1 Eq 49
res_mat_2=zeros(num1_to_search,num2_to_search);%diff2 Eq 50
res_mat_3=zeros(num1_to_search,num2_to_search);%abs(diff1*diff2)

tol=0.05;
for i=1:num1_to_search
    for j=1:num2_to_search
       [tmp1,tmp2]=SolTwoEquOU( Q_h_vector(i),Q_l_vector(j));
%  [tmp1,tmp2]=SolTwoEqu( Q_h_vector(i),Q_l_vector(j));
        res_mat_1(i,j)=tmp1;
        res_mat_2(i,j)=tmp2;
        res_mat_3(i,j)=abs(tmp1)+abs(tmp2);
%         if res_mat_1(i,j)<tol && res_mat_1(i,j)>-tol  && res_mat_1(i,j)<tol && res_mat_1(i,j)>-tol 
%             Q_h_vector(i),Q_l_vector(j)
%             break
%         end
    end
end
min_val_1=min(min(res_mat_3))    
[row,col] = find(res_mat_3==min(res_mat_3(:)));
Q_h_sol=Q_h_vector(row)
Q_l_sol=Q_l_vector(col)

Diff1_sol=res_mat_1(row,col)
Diff2_sol=res_mat_2(row,col)
