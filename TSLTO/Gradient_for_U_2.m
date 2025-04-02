function [value,a] = Gradient_for_U_2(G,U1,U2,U3,L,beta,rou,Y,T,V,output_mode)
Delta = G;
% Delta = nmodeproduct(Delta,U1,1);
% Delta = nmodeproduct(Delta,U2,2);
% Delta = nmodeproduct(Delta,U3,3);
% Delta = Delta - L;
Delta = ttm(Delta,{U1,U2,U3},[1,2,3])-L;
v = unfold(beta.*Delta,2)*kron(U3,U1)*unfold(G,2)';
%v = ten2mat(beta.*Delta,size(Delta),2)*kron(U3,U1)*ten2mat(G,size(G),2)';
v = v - T{2}'*V{2};
%v = (v + rou.*T{2}'*(T{2}*U2-Y{2}));
v = (v + rou.*T{2}'*(diff(U2,1,1)-Y{2}));
%value = k1+k2-k3+k4;
%value = v - U2*Sym(U2'*v);
a = v*U2' - U2*v';
value = a*U2;
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
