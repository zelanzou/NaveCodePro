function [gamma,Wm,Wc] = UKFParameter2(n)
beta = 2; kapa = 0;alpha = 0.001;
lamda = alpha^2*(n+kapa)-n;
gamma = sqrt( n+ lamda);
m = gamma^2;
Wm = [lamda/m; repmat(1 /(2*m),2*n,1)];
Wc = [lamda/m + 1 +beta - alpha^2;Wm(2:end)];
end