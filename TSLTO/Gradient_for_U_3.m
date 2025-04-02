function [value,a] = Gradient_for_U_3(G,U1,U2,U3,L,beta,rou,Y,T,V,output_mode)
Delta = G;
% Delta = nmodeproduct(Delta,U1,1);
% Delta = nmodeproduct(Delta,U2,2);
% Delta = nmodeproduct(Delta,U3,3);
% Delta = Delta - L;
Delta = ttm(Delta,{U1,U2,U3},[1,2,3])-L;
v = unfold(beta.*Delta,3)*kron(U2,U1)*unfold(G,3)';
%v = ten2mat(beta.*Delta,size(Delta),3)*kron(U2,U1)*ten2mat(G,size(G),3)';
v = v - T{3}'*V{3};
%v = (v + rou.*T{3}'*(T{3}*U3-Y{3}));
v = (v + rou.*T{3}'*(diff(U3,1,1)-Y{3}));
%value = k1+k2-k3+k4;
%value = v - U3*Sym(U3'*v);
a = v*U3' - U3*v';
value = a*U3;
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
