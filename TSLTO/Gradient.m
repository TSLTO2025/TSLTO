function grad = Gradient(R,T,T_r,W,Z,gamma,X,L,k,P)
%M1^(-1)(Tr'*W*Tl)
Tl_W_Tr = T{1}'*W*T_r;
%Tl_W_Tr = diff(diff(W,1,1),1,2);
sizeR = size(R);
%M1^(-1):将矩阵沿mode-1折回到张量形式
%grad = - mat2ten(Tl_W_Tr,sizeR,1);
grad = - fold(Tl_W_Tr,1,sizeR);
%gamma*M1^(-1)(Tl*R*Tr')
%Tl_R_Tr = T{1}*ten2mat(R,sizeR,1)*T_r';
%Tl_R_Tr = T{1}*unfold(R,1)*T_r';
Tl_R_Tr = diff(diff(unfold(R,1),1,1),1,2);
%grad = grad + gamma.*mat2ten(T{1}'*(Tl_R_Tr - Z)*T_r,sizeR,1);
grad = grad - gamma.*fold(T{1}'*(Z - Tl_R_Tr)*T_r,1,sizeR);
%grad = grad + gamma.*fold(diff(diff((Tl_R_Tr - Z),1,1),1,2),1,sizeR);
grad = grad - P + k.*(R+L-X);
end