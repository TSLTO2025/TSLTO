function F = curr_value1(G,U,L,R,T,T_r,lambda,mu,beta)
approx_tensor = double(ttm(G,U,[1,2,3])) - L;
F = norm((approx_tensor),'fro')^2;
end