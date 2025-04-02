function F3 = curr_value2(G,U,L,R,T,T_r,lambda,mu,beta)
F3 = 0;%l2,0范数部分
for i = 1:3
    Ti_Ui = T{i}*U{i};
    F3 = F3 + nnz(round(sqrt(sum((Ti_Ui)'.^2,1)),3));
end