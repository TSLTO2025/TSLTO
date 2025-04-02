function T = Toeplitz(n)
%input: rows
%output:toeplitz matrix \in R^{(n-1) \times n}
%initiallize
T = zeros(n-1,n);
for i = 1:(n-1)
    T(i,i) = 1;
    T(i,i+1) = -1;
end