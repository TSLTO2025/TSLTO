function [X,L,R,U,iter,Y,G] = main(RR,RealT,Tensor,Omega,tucker_size,alpha,beta,gamma,lambda,mu,maxiter,maxiter0,epsilon,epsilon_0,delta,ita,ro,tau1,tau2,lambda_0,rho,k)
%input
%Tensor:原观测张量 %Omega: %tucker_size:e.g.[4 5 6] %alpha:3*1 %lambda:3*1 %mu:2*1 %maxiter:最大迭代次数
%output
%X:补全后的张量 %W:低秩张量 %R:异常值

% t1 = zeros(maxiter,1);
% t2 = zeros(maxiter,1);
% x = 3:maxiter;
% y1 = zeros(1, maxiter-2); 
% y2 = zeros(1, maxiter-2); 
% y3 = zeros(1, maxiter-2); 
% figure;
% subplot(3, 1, 1);
% h1 = plot(x, y1, 'r');
% title('norm L');
% xlabel('迭代次数');
% ylabel('值');
% 
% subplot(3, 1, 2);
% h2 = plot(x, y2, 'g');
% title('norm R');
% xlabel('迭代次数');
% ylabel('值');
% 
% subplot(3, 1, 3);
% h3 = plot(x, y3, 'b');
% title('norm X');
% xlabel('迭代次数');
% ylabel('值');

%初始化
sizeT = size(Tensor);
%anomaly tensor R
R = zeros(sizeT);
%R = Tensor - reconstructTensor(G,U);
%R = double(Tensor - ttm(G,U,[1,2,3]));
%R = RR;
% R(~Omega) = -mean(Tensor(:));
%tucker_size = floor(tucker_size .* sizeT);
tucker_Tensor = tucker_als(tensor(Tensor),tucker_size);
%core tensor G
G = tucker_Tensor.core;
%factor matrices U
U = tucker_Tensor.U;
for i = 1:3
    %U{i} = round(U{i},3);
    U{i} = toep_sparse_column_ortho(sizeT(i),tucker_size(i),5);
end
%V
V = cell(3,1);
for i = 1:3
    V{i} = rand(sizeT(i)-1,tucker_size(i));
end
%T 1,2,3
T = cell(3,1);
for i = 1:3
    T{i} = sparse(Toeplitz(sizeT(i)));
end
%T_l = T{1};
T_r = sparse(Toeplitz(prod(sizeT)/sizeT(1)));
%Y
Y = cell(3,1);
for i = 1:3
    %Y{i} = T{i}*U{i};
    Y{i} = rand(sizeT(i)-1,tucker_size(i));
end
%Z
Z = rand([sizeT(1)-1 prod(sizeT)/sizeT(1)-1]);
%Z = T{1}*unfold(R,1)*T_r';
%W
W = rand(sizeT(1)-1,prod(sizeT)/sizeT(1)-1);
%P
P = rand(sizeT);
%L = rand(sizeT);
L = Tensor;
%Omega_c = 1-Omega;
%prev_f = curr_value(G,U,L,R,T,T_r,lambda,mu,beta);
X_pre = Tensor;
G_pre = G;
% U_pre = U;
R_pre = R;
L_pre = L;
% Y_pre = Y;
% Z_pre = Z;
% V_pre = V;
% W_pre = W;
% P_pre = P;
% value1 = zeros(maxiter,1);
% value2 = zeros(maxiter,1);
% value3 = zeros(maxiter,1);
% value4 = zeros(maxiter,1);
% value5 = zeros(maxiter,1);

for iter = 1:maxiter
    %更新X
    X = Tensor.*Omega + (L+R-(1/k).*P).*(1 - Omega);
    dif_X = norm(X_pre - X,'fro')/norm(X_pre,'fro');

    %更新G
    G = L;
    % for i = 1:3
    %     G = nmodeproduct(G,(U{i})',i);
    % end
    G = double(ttm(tensor(G),{U{1}',U{2}',U{3}'},[1,2,3]));
    dif_G = norm(double(G_pre-G),'fro')/norm(double(G_pre),'fro');

    %更新Yi
    for i = 1:3
        b = sqrt(2*lambda(i)/alpha(i));
        %x = T{i}*U{i} - V{i}./alpha(i);
        x = diff(U{i},1,1) - (1/alpha(i)).*V{i};
        aa = zeros(size(V{i}));
        for j = 1:(sizeT(i)-1)
            if norm(x(j,:)) > b
                aa(j,:) = x(j,:);
            end
        end
        Y{i} = aa;
    end
    % dif_Y1 = norm(Y_pre{1} - Y{1},'fro')/norm(Y_pre{1},'fro');
    % dif_Y2 = norm(Y_pre{2} - Y{2},'fro')/norm(Y_pre{2},'fro');
    % dif_Y3 = norm(Y_pre{3} - Y{3},'fro')/norm(Y_pre{3},'fro');

    %更新U
    maxiter_1 = 100;
    maxiter_2 = 50;
    U{1} = Curv_1(tensor(G),U{1},U{2},U{3},L,beta,alpha(1),Y,T,V,tau1,tau2,delta,ita,epsilon_0,maxiter_1,maxiter_2,ro);
    U{2} = Curv_2(tensor(G),U{1},U{2},U{3},L,beta,alpha(2),Y,T,V,tau1,tau2,delta,ita,epsilon_0,maxiter_1,maxiter_2,ro);
    U{3} = Curv_3(tensor(G),U{1},U{2},U{3},L,beta,alpha(3),Y,T,V,tau1,tau2,delta,ita,epsilon_0,maxiter_1,maxiter_2,ro);
    % dif_U1 = norm(U_pre{1} - U{1},'fro')/norm(U_pre{1},'fro');
    % dif_U2 = norm(U_pre{2} - U{2},'fro')/norm(U_pre{2},'fro');
    % dif_U3 = norm(U_pre{3} - U{3},'fro')/norm(U_pre{3},'fro');

    %更新R
    %R = PGM(R, T, T_r, W, Z, X, L, gamma, k, P,mu,lambda_0,rho,maxiter);
    [R,lambda_0] = PGM(R, T, T_r, W, Z, X, L, gamma, k, P,mu,lambda_0,rho,20);
    %R = PGM(R, T, T_r, W, Z, X, L, gamma, k, P,mu,rho,100);
    dif_R = norm(double(R_pre - R),'fro')/norm(double(R_pre),'fro');
    if iter == 25000
        lambda_0 = lambda_0 * 0.8;
    end
    if iter == 30000
        lambda_0 = lambda_0 * 0.6;
    end
    if iter == 31000
        lambda_0 = lambda_0 * 0.5;
    end
    if iter == 35000
        lambda_0 = lambda_0 * 0.5;
    end
    if iter == 40000
        lambda_0 = lambda_0 * 0.5;
    end
    if iter == 50000
        lambda_0 = lambda_0 * 0.4;
    end
    if iter == 55000
        lambda_0 = lambda_0 * 0.4;
    end

    %更新L
    Delta = G;
    % Delta = nmodeproduct(Delta,U{1},1);
    % Delta = nmodeproduct(Delta,U{2},2);
    % Delta = nmodeproduct(Delta,U{3},3);
    Delta = double(ttm(tensor(Delta),U,[1,2,3]));
    L = (beta/(beta+k)).*Delta + (k/(beta+k)).*(X-R) + (1/(beta+k)).*P;
    dif_L = norm(double(L_pre - L),'fro')/norm(double(L_pre),'fro');

    % if iter >=3
    %     y2(iter-2) = norm(R,'fro');
    %     subplot(3, 1, 2);
    %     set(h1, 'YData', y1);
    % 
    %     y1(iter-2) = norm(L,'fro');
    %     subplot(3, 1, 1);
    %     set(h2, 'YData', y2);
    % 
    %     y3(iter-2) = norm(X,'fro');
    %     subplot(3, 1, 3);
    %     set(h3, 'YData', y3);
    %     pause(0.1);
    % end

    %更新Z
    c = 2*mu(2)/gamma;
    % x = T{1}*ten2mat(R,sizeT,1)*T_r' - (1/gamma) .* W;
    % x = T{1}*unfold(R,1)*T_r' - (1/gamma) .* W;
    x = diff(diff(unfold(R,1),1,1),1,2) - (1/gamma) .* W;
    %x = diff(diff(ten2mat(R,sizeT,1),1,1),1,2) - (1/gamma) .* W;
    Z = prox_zero(c,x);
    % Z = zeros(sizeT(1)-1,prod(sizeT)/sizeT(1)-1);
    % for i = 1:size(W,1)
    %     for j = 1:size(W,2)
    %         if x(i,j)^2 > c
    %             Z(i,j) = x(i,j);
    %         end
    %     end
    % end
    % dif_Z = norm(Z_pre - Z,'fro')/norm(Z_pre,'fro');

    %更新Vi
    for i = 1:3
        %V{i} = V{i} + alpha(i).*(Y{i} - T{i}*U{i});
        V{i} = V{i} + alpha(i).*(Y{i} - diff(U{i},1,1));
    end
    % dif_V1 = norm(V_pre{1} - V{1},'fro')/norm(V_pre{1},'fro');
    % dif_V2 = norm(V_pre{2} - V{2},'fro')/norm(V_pre{2},'fro');
    % dif_V3 = norm(V_pre{3} - V{3},'fro')/norm(V_pre{3},'fro');

    %更新W
    %W = W + gamma.*(Z - T{1}*ten2mat(R,sizeT,1)*T_r');
    % W = W + gamma.*(Z - T{1}*unfold(R,1)*T_r');
    W = W + gamma.*(Z - diff(diff(unfold(R,1),1,1),1,2));
    %W = W + gamma.*(Z - diff(diff(ten2mat(R,sizeT,1),1,1),1,2));
    % dif_W = norm(W_pre - W,'fro')/norm(W_pre,'fro');

    %更新P
    P = P + k.*(X - L - R);
    %dif_P = norm(double(P_pre - P),'fro')/norm(double(P_pre),'fro');

    %收敛条件
    %curr_f = curr_value(G,U,L,R,T,T_r,lambda,mu,beta);
    %diff = abs(curr_f-prev_f)/(abs(prev_f)+1);
    %t1(iter) = diff;
    %dif = norm(X_pre - X,'fro')/norm(X_pre,'fro');
    %t2(iter) = dif;
    % if iter>1 && (diff < epsilon) && (dif < epsilon)
    %     disp('已达收敛条件');
    %     break;
    % end
    % fprintf('difX = %d', dif_X);
    % fprintf('difG = %d', dif_G);
    % fprintf('difR = %d', dif_R);
    % fprintf('difL = %d', dif_L);
    if iter>1 && (dif_X < epsilon) && (dif_G < epsilon)&& (dif_R < epsilon) && (dif_L < epsilon)
        disp('已达收敛条件');
        break;
    end
    X_pre = X;
    G_pre = G;
    %U_pre = U;
    R_pre = R;
    L_pre = L;
    % Y_pre = Y;
    % Z_pre = Z;
    % V_pre = V;
    % W_pre = W;
    % P_pre = P;
    %fprintf('momo: iter = %d\n   f value=%f\n', iter, curr_f);
    %prev_f = curr_f;

    %fprintf('main: iterations = %d', iter);
    if mod(iter, 10) == 0
        fprintf('\ndifX = %d\n', dif_X);
        fprintf('difG = %d\n', dif_G);
        fprintf('difR = %d\n', dif_R);
        fprintf('difL = %d\n', dif_L);
        fprintf('main: iterations = %d\n', iter);
        l = double(ttm(tensor(G),U,[1,2,3]));
        tensor_errors(RealT,l,Omega);
        tensor_errors(RealT,L,Omega);
        %tensor_errors(RealT,X,Omega);
        anomaly_test(RR,R);
    end
    % l = double(ttm(tensor(G),U,[1,2,3]));
    % tensor_errors(RealT,l,Omega);
    % tensor_errors(RealT,L,Omega);
    % anomaly_test(RR,R);

    %更新alpha,gamma
    alpha = min(alpha.*1.15,1e4);
    gamma = min(gamma*1.15,1e4);
    %beta = min(beta*1.15,10);
    k = min(k*1.15,1e4);
    %mu(1) = max(mu(1)*0.8,1e-5);
    % if iter <= 15
    %     mu(1) = max(mu(1)*0.8,1e-5);
    % end
    %mu(2)=min(mu(2)*1.15,1e3);
    % if iter > 30
    %     mu(1) = min(mu(1)*1.05,10);
    % end
    G = tensor(G);
    % value1(iter) = curr_value1(G,U,L,R,T,T_r,lambda,mu,beta);
    % value2(iter) = curr_value2(G,U,L,R,T,T_r,lambda,mu,beta);
    % value3(iter) = curr_value3(G,U,L,R,T,T_r,lambda,mu,beta);
    % value4(iter) = curr_value4(G,U,L,R,T,T_r,lambda,mu,beta);
    % value5(iter) = curr_value(G,U,L,R,T,T_r,lambda,mu,beta);
end
fprintf('main ends: total iterations = %d\n', iter);
% val = cell(6,1);
% val{1} = value1(1:iter);
% val{2} = value2(1:iter);
% val{3} = value3(1:iter);
% val{4} = value4(1:iter);
% val{5} = value5(1:iter);
