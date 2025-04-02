function tau = tau_1(s,y)
tau = sum(dot(s,s))/abs(sum(dot(s,y)));
end