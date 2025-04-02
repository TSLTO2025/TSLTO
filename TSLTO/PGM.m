%function x_0 = PGM(x_0, T, T_r, W, Z, X, L, gamma, k, P,mu,lambda_0,rho,maxiter)
function [x_0,lambda] = PGM(x_0, T, T_r, W, Z, X, L, gamma, k, P,mu,lambda_0,rho,maxiter)
lambda = lambda_0;
%x = x_0;
g = Gradient(x_0,T,T_r,W,Z,gamma,X,L,k,P);
for i=1:maxiter
    a = 2*lambda*mu(1);
    y = double(x_0-lambda.*g);
    b = prox_zero(a,y);
    kk = original_f_value(b, T, T_r, W, Z, X, L, gamma, k, P);
    % if mod(i, 1) == 0
    %     %fprintf('PGM: x_0-lambda.*g = %d',y);
    %     fprintf('PGM: iterations = %d value=%f\n', i, kk);
    % end
    if kk<=f_hat(x_0, T, T_r, W, Z, X, L, gamma, k, P, b, lambda)
        x_0 = b;
        break;
    else
    lambda = rho*lambda;
    end
end
