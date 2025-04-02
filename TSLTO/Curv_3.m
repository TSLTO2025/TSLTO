function X_pre = Curv_3(G,U1,U2,X_0,L,beta,alpha,Y,T,V,tau1,tau2,delta,ita,epsilon,maxiter_1,maxiter_2,ro)
%X_0: 初始点 即U1
%tau1: min tau
%tau2: max tau
%maxiter_1: for inner loop
%maxiter_2: for outer loop
t_1 = zeros(maxiter_2,1);
t_2 = zeros(maxiter_2,1);
X_pre = X_0;
q = 1;
c = f_U3(G,U1,U2,X_pre,L,beta,alpha,Y,T,V);
tau = tau1;
F_pre = c;
for i = 1:maxiter_2
    %X_pre = X_cur;
    [g_pre,a] = Gradient_for_U_3(G,U1,U2,X_pre,L,beta,alpha,Y,T,V,3);
    if (norm(g_pre,'fro') <= epsilon) && (i >=5) 
        %disp(i);
        break;
    end
    for j = 1:maxiter_1
        tau = delta*tau;
        X_cur = Retr_Y(X_pre,a,tau);
        F = f_U3(G,U1,U2,X_cur,L,beta,alpha,Y,T,V);
        if F < c - 1/2*ro*tau*norm(a,"fro")^2
            break;
        end
    end
    g_cur = Gradient_for_U_3(G,U1,U2,X_cur,L,beta,alpha,Y,T,V,1);
    q = ita*q + 1;
    c = (1/q).*(ita*q*c + F);
    s = X_cur - X_pre;
    tol_1 = norm(s,'fro')/sqrt(size(X_0,1));
    t_1(i) = tol_1;
    tol_2 = (F_pre - F)/(1 + abs(F_pre));
    t_2(i) = tol_2;
    if (tol_1 < 1e-3) && (tol_2 < 1e-6) && (i >=5)
        % disp(i);
        % disp(2);
        break;
    end
    if (i > 10) && (mean(t_1(i-9:i)) <= 1e-2) && (mean(t_2(i-9:i)) <= 1e-5)
        % disp(i);
        % disp(3);
        break;
    end
    y = g_cur - g_pre;
    tau = max(tau1,min(tau_1(s,y),tau2));
    X_pre = X_cur;
    % if mod(i, 1e4) == 0  %每过20次iteration，就print一次结果
    %     fprintf('Curv_3: iterations = %d   value=%f\n', i, F_pre);
    % end
    F_pre = F;
end
%fprintf("fU3 final: %f\n  iter: %d\n",F_pre,i);