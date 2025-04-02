function F = curr_value(G,U,L,R,T,T_r,lambda,mu,beta)
approx_tensor = double(ttm(G,U,[1,2,3])) - L;
F = 0.5*beta*norm((approx_tensor),'fro')^2;
F3 = 0;%l2,0范数部分
for i = 1:3
    %Ti_Ui = T{i}*U{i};
    Ti_Ui = diff(U{i},1,1);
    F3 = F3 + lambda(i)*nnz(round(sqrt(sum((Ti_Ui)'.^2,1)),3));
end
F = F+F3;
F = F+mu(1)*nnz(round(R,6)); 
%TRT = T{1}*ten2mat(R,size(R),1)*T_r'; 
%TRT = T{1}*unfold(R,1)*T_r'; 
TRT = diff(diff(unfold(R,1),1,1),1,2);
%TRT = diff(diff(ten2mat(R,size(R),1),1,1),1,2);
F = F+mu(2)*nnz(round(TRT,6)); 
end