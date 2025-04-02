function ita = Backtracking_armijo_2(x,t_0,c,tau,maxiter,G,Omega,alpha,U1,U3,R,Tensor,beta,Y,T,V)
t = t_0;
ita_0 = -Gradient_for_U_2(G,Omega,U1,x,U3,R,Tensor,beta,alpha,Y,T,V,1);
for i = 1:maxiter
    t = t*tau;
    
    a = f_U2(G,Omega,U1,x,U3,R,Tensor,beta,alpha,Y,T,V) - f_U2(G,Omega,U1,Retr_polar(x,t.*ita_0),U3,R,Tensor,beta,alpha,Y,T,V);
    b = c*t*sum(dot(-Gradient_for_U_2(G,Omega,U1,x,U3,R,Tensor,beta,alpha,Y,T,V,1),ita_0));
    if a >= b
        break;
    end
end
ita = t.*ita_0;