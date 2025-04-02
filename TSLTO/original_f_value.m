function f_Z = original_f_value(R, T, T_r, W, Z, X, L, gamma, k, P)
% Z_TRT = Z-T{1}*ten2mat(R,size(R),1)*T_r';
%Z_TRT = Z-T{1}*unfold(R,1)*T_r';
Z_TRT = Z-diff(diff(unfold(R,1),1,1),1,2);
%Z_TRT = Z-diff(diff(ten2mat(R,size(R),1),1),1,2);
f_Z = dot(Z_TRT(:),W(:));
f_Z = f_Z + (gamma/2)*norm(Z_TRT,'fro')^2;
X_L_R = double(X-L-R);
f_Z = f_Z + dot(X_L_R(:),P(:));
f_Z = f_Z + (0.5*k)*norm(X_L_R,'fro')^2;
end