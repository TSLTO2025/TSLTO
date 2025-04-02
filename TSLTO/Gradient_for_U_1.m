function [value,a] = Gradient_for_U_1(G,U1,U2,U3,L,beta,alpha,Y,T,V,output_mode)
Delta = tensor(G);
% Delta = nmodeproduct(Delta,U1,1);
% Delta = nmodeproduct(Delta,U2,2);
% Delta = nmodeproduct(Delta,U3,3);
% Delta = Delta - L;
Delta = ttm(Delta,{U1,U2,U3},[1,2,3])-L;%ttm代替nmodeproduct，更快
v = unfold(beta.*Delta,1)*kron(U3,U2)*unfold(G,1)';
%v = ten2mat(beta.*Delta,size(Delta),1)*kron(U3,U2)*ten2mat(G,size(G),1)';
v = v - T{1}'*V{1};
%v = (v + alpha.*T{1}'*(T{1}*U1-Y{1}));
v = (v + alpha.*T{1}'*(diff(U1,1,1)-Y{1}));
%value = k1+k2-k3+k4;
%value = v - U1*Sym(U1'*v);
a = v*U1' - U1*v';
value = a*U1;
% 根据 output_mode 决定输出
if output_mode == 1
    % 只输出 value
    a = []; % 不输出 a
elseif output_mode == 2
    % 只输出 a
    value = []; % 不输出 value
elseif output_mode == 3
    % 同时输出 value 和 a
    % 不需要做任何修改，直接返回 value 和 a
else
    error('Invalid output_mode. Use 1 for value only, 2 for a only, or 3 for both.');
end
end