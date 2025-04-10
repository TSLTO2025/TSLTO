%GUANGZHOU 
Tensor = tensor;
clear tensor;
Tensor = permute(Tensor,[1,3,2]);
%%
hosvd(Tensor,0.008);
tt = tucker_als(Tensor,[2,5,6]);
tensor_errors(double(Tensor),double(ttm(tt.core,tt.U,[1,2,3])),zeros(size(Tensor)));
%%
rng(42);
tuckerT1 = tucker_als(tensor(Tensor-20),[2,5,6]);
U1 = tuckerT1.U;
G1 = tuckerT1.core;
g = G1.*1.4e-1;
t = double(ttm(g,U1,[1,2,3]));
% t = double(Tensor-20).*0.12;
% figure;
% imagesc(Tensor(:,:,6));
% colorbar;
%%
rng(42);
O = zeros(214, 8784);
block_size = [150,25]; % 每个块的大小
a = randi(64,50,1);
b = randi(8759,50,1);
% block_size = [3,783]; % 每个块的大小
% a = randi(211,80,1);
% b = randi(8001,80,1);
block_positions = [a,b];
% 将每个块填充为1
for i = 1:size(block_positions, 1)
    row_start = block_positions(i, 1);
    col_start = block_positions(i, 2);
    row_end = row_start + block_size(1) - 1;
    col_end = col_start + block_size(2) - 1;
    O(row_start:row_end, col_start:col_end) = 1;
end
%
% R = zeros(50,2500);
% R(O == 1) = randn(nnz(O), 1).*2;
R = zeros(214,144*61);
R(O == 1) = (0.1.*randn(nnz(O), 1)+ 4);
figure;
imagesc(R);
colorbar;
RR = fold(R,1,size(Tensor));
%%
% Tensor = permute(Tensor,[1,3,2]);
%Tensor = Tensor(1:30,:,1:30);
%RR = RR(1:30,:,1:30);
B = t + RR;
figure;
imagesc(B(:,:,18));
colorbar;
%% 
rng(42);
[Observed_Tensor,Omega] = random_missing(B,0.3);
tucker_size = [2 5 6]; 
alpha = [2e1,2e1,2e1];
k = 0.3;
lambda = [1.2,1.2,1.2]; 
mu = [5,20]; 
beta = 450; 
gamma = 10; 
maxiter = 200000; 
maxiter0 = 60;
epsilon = 1e-5;
epsilon_0 = 1e-4;
rng(42);
tic;
[X,L,R,U,iter,Y,G] = main(RR,t,Observed_Tensor,Omega,tucker_size,alpha,beta,gamma,lambda,mu,maxiter,maxiter0,epsilon,epsilon_0,0.3,0.01,0.3,1e-2,0.9,3,0.8,k);
running_time = toc;
%%
tensor_errors(B,X,Omega);
tensor_errors(t,L,Omega);
tensor_errors(t,double(ttm(tensor(G),U,[1,2,3])),Omega);
anomaly_test(RR,R);
fprintf('还原后：\n');
abc = double(ttm(tensor(G./0.14),U,[1,2,3]));
smoothT = double(ttm(G1,U1,[1,2,3]));
tensor_errors(smoothT+20,abc+20,Omega);
%%
slice = 8; 
figure;
subplot(1,8,1);
imagesc(t(:,:,slice));
title('原始低秩张量切片');
colorbar;
subplot(1,8,2);
imagesc(RR(:,:,slice));
title('原始稀疏张量切片');
colorbar;
subplot(1,8,3);
imagesc(B(:,:,slice));
title('原始张量切片');
colorbar;
subplot(1,8,4);
imagesc(Observed_Tensor(:,:,slice));
title('观测张量切片');
colorbar;
tucker_part = double(ttm(tensor(G),U,[1,2,3]));
X = double(X);
L = double(L);
R = double(R);
subplot(1,8,5);
imagesc(X(:,:,slice));
title('补全后张量切片X');
colorbar;
subplot(1,8,6);
imagesc(L(:,:,slice));
title('低秩张量切片L');
colorbar;
subplot(1,8,7);
imagesc(R(:,:,slice));
title('块稀疏张量切片R');
colorbar;
subplot(1,8,8);
imagesc(tucker_part(:,:,slice));
title('tucker低秩部分');
colorbar;
