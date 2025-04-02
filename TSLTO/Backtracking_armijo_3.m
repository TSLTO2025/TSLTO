function ita = Backtracking_armijo_3(x,t_0,c,tau,maxiter,G,Omega,alpha,U1,U2,R,Tensor,beta,Y,T,V)
t = t_0;
ita_0 = -Gradient_for_U_3(G,Omega,U1,U2,x,R,Tensor,beta,alpha,Y,T,V,1);
for i = 1:maxiter
    t = t*tau;
    
    a = f_U3(G,Omega,U1,U2,x,R,Tensor,beta,alpha,Y,T,V) - f_U3(G,Omega,U1,U2,Retr_polar(x,t.*ita_0),R,Tensor,beta,alpha,Y,T,V);
    b = c*t*sum(dot(-Gradient_for_U_3(G,Omega,U1,U2,x,R,Tensor,beta,alpha,Y,T,V,1),ita_0));
    if a >= b
        break;
    end
end
ita = t.*ita_0;