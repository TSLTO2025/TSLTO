%模拟数据低秩部分
rng(42);
r = 3;
G = tensor(rand(r,r,r)).*100;
%G = tensor(rand(r,r,r));
u = cell(3,1);
for i = 1:3
    %U{i} = orth(randn(100,r))*diag([1:r].^(-0.5));
    u{i} = toep_sparse_column_ortho(50,3,2);
end
Tensor = double(ttm(G,u,[1,2,3]));
figure;
imagesc(Tensor(:,:,12));
colorbar;
%%
%模拟数据稀疏部分
rng(42);
O = zeros(50, 2500);
block_size = [2,125]; %每个块的大小
%手动指定15个块的位置
a = randi(48,50,1);
b = randi(2375,50,1);
block_positions = [a,b];
%将每个块填充为1
for i = 1:size(block_positions, 1)
    row_start = block_positions(i, 1);
    col_start = block_positions(i, 2);
    row_end = row_start + block_size(1) - 1;
    col_end = col_start + block_size(2) - 1;
    O(row_start:row_end, col_start:col_end) = 1;
end
R = zeros(50,2500);
R(O == 1) = (0.1.*randn(nnz(O), 1)+ 4); %均值为4，方差为1
% R = ones(size(O)).*4; % 生成与 O 大小相同的随机矩阵
% R(~O) = 0;
figure;
imagesc(R);
colorbar;
RR = fold(R,1,[50,50,50]);
%%
rng(42);
B=Tensor+RR;
[Observed_Tensor,Omega] = random_missing(B,0.3);
tucker_size = [3 3 3];
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
%%
rng(42);
tic;
[X,L,R,U,iter,Y,G] = main(RR,Tensor,Observed_Tensor,Omega,tucker_size,alpha,beta,gamma,lambda,mu,maxiter,maxiter0,epsilon,epsilon_0,0.3,0.01,0.3,1e-2,0.9,3,0.8,k);
running_time = toc;
tucker_part = double(ttm(tensor(G),U,[1,2,3]));
tensor_errors(Tensor,tucker_part,Omega);
tensor_errors(Tensor,L,Omega);
anomaly_test(RR,R);
fprintf('time: %f\n',running_time);
%%
slice = 12; %一个切片可视化
figure;
subplot(1,8,1);
imagesc(Tensor(:,:,slice));
title('原始低秩张量切片');
colorbar;
subplot(1,8,2);
imagesc(RR(:,:,slice));
title('原始稀疏张量切片');
colorbar;
subplot(1,8,3);
imagesc(Observed_Tensor(:,:,slice));
title('观测张量切片');
colorbar;
subplot(1,8,4);
imagesc(X(:,:,slice));
title('补全后张量切片X');
colorbar;
subplot(1,8,5);
imagesc(L(:,:,slice));
title('低秩张量切片L');
colorbar;
subplot(1,8,6);
r = L - tucker_part;
r = round(r,3);
imagesc(R(:,:,slice));
title('块稀疏张量切片R');
colorbar;
subplot(1,8,7);
imagesc(tucker_part(:,:,slice));
title('tucker分解得到的低秩张量切片');
colorbar;
subplot(1,8,8);
imagesc(r(:,:,slice));
title('L - tucker');
colorbar;
%%
B = Tensor + RR;
figure('Color','white');
subplot(1,3,1);
imagesc(B(:,:,1));
title('(a)');
colorbar;
subplot(1,3,2);
imagesc(Tensor(:,:,1));
title('(b)');
colorbar;
subplot(1,3,3);
imagesc(RR(:,:,1));
title('(c)');
colorbar;