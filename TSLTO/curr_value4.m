function F = curr_value4(G,U,L,R,T,T_r,lambda,mu,beta)
%TRT = T{1}*ten2mat(R,size(R),1)*T_r'; 
%TRT = T{1}*unfold(R,1)*T_r'; 
TRT = diff(diff(unfold(R,1),1,1),1,2);
%TRT = diff(diff(ten2mat(R,size(R),1),1,1),1,2);
F = nnz(round(TRT,6)); 
end