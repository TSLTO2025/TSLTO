function Q = toep_sparse_column_ortho(K,I,s_row)
%K:行数
%I:列数
%s_row:非toep-l20稀疏行数
U = randn(K,I);
random_numbers = randperm(K,K-s_row);
%generate zero rows
for i = 2:(K-s_row)
    U(random_numbers(i),:) = U(random_numbers(1),:);
end
%Q,R decomposition
[Q,~] = qr(U);
Q     = Q(:,1:I);
end