function x = prox_zero(k,y)
%k:参数
%Y:括号内的值(matrix)
%k:参数
%Y:括号内的值(matrix)
x = zeros(size(y));
x(y.*y>k) = y(y.*y>k);
