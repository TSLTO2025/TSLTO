function y = Retr_Y(X,A,tau)
I = eye(size(A,1));
% %求W
% max_a = -1;
% ii = 0;
% for i = 1:size(X,2)
%     a = trace(X*G'*(G(:,i)*X(:,i)' - X(:,i)*G(:,i)'));
%     if a > max_a
%         ii = i;
%     end
% end
% W = G(:,ii)*X(:,ii)' - X(:,ii)*G(:,ii)';
%最终结果
y = (I + tau/2.*A)\((I-tau/2.*A)*X);
