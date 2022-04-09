function [gamma,Wm,Wc] = UKFParameter(n)
beta = 2; kapa = 0;alpha = 0.001;
lamda = alpha^2*(n+kapa)-n;
gamma = sqrt( n+ lamda);
m = gamma^2;
Wm = [lamda/m; repmat(1 /(2*m),2*n,1)];
Wc = eye(2*n+1)*1 /(2*m);
Wc(1,1) = lamda/m + 1 +beta - alpha^2;
end