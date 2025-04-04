function value = f_U3(G,U1,U2,U3,L,beta,alpha,Y,T,V)
Delta = tensor(G);
% Delta = nmodeproduct(Delta,U1,1);
% Delta = nmodeproduct(Delta,U2,2);
% Delta = nmodeproduct(Delta,U3,3);
% Delta = Delta - L;
Delta = double(ttm(Delta,{U1,U2,U3},[1,2,3]))-L;
value = (beta/2) * (norm(Delta,'fro'))^2;
t= Y{3}-diff(U3,1,1);
value = value + sum(dot(t,V{3}));
value = value + alpha/2 * norm(t,"fro")^2;
