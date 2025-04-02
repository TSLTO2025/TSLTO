function ita = Backtracking_armijo_1(x,t_0,c,tau,maxiter,G,Omega,alpha,U2,U3,R,Tensor,beta,Y,T,V)
t = t_0;
ita_0 = -Gradient_for_U_1(G,Omega,x,U2,U3,R,Tensor,beta,alpha,Y,T,V,1);
for i = 1:maxiter
    t = t*tau;
    a = f_U1(G,Omega,x,U2,U3,R,Tensor,beta,alpha,Y,T,V) - f_U1(G,Omega,Retr_polar(x,t.*ita_0),U2,U3,R,Tensor,beta,alpha,Y,T,V);
    b = c*t*sum(dot(-Gradient_for_U_1(G,Omega,x,U2,U3,R,Tensor,beta,alpha,Y,T,V,1),ita_0));
    if a >= b
        break;
    end
end
ita = t.*ita_0;