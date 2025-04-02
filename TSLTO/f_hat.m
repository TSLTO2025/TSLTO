function f_h = f_hat(x_0, T, T_r, W, Z, X, L, gamma, k, P, b, lambda)
%f_hat(b,x)=f(x)+(grad_f(x)^T)*(b-x)+(1/(2*λ))*||b-x||_F^2
%f(x)
f_h = original_f_value(x_0, T, T_r, W, Z, X, L, gamma, k, P);
%grad_f(x)
grad_f_x = Gradient(x_0,T,T_r,W,Z,gamma,X,L,k,P);
%(grad_f(x)^T)*(b - x)
dif = double(b-x_0);
f_h = f_h + dot(grad_f_x(:), dif(:));
%(1/(2*λ))*||b-x||_2^2
f_h = f_h + (1/(2*lambda))*(norm(dif,'fro'))^2;
end