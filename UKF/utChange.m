function [X0, P0] = utChange(xk,Pk,Wm,Wc,tpara)
n = size(kf.Pk,1);
L = chol( Pk,'lower');
x = repmat(xk,1,n);
xsigma0 = [xk,x+gamma*L,x-gamma*L];%15*31
Y(:,1)  = feval(hfx, xsigma0(:,1), tpara);
m=length(Y); y = Wm(1)*Y(:,1);
Y = repmat(Y,1,2*n+1);%∑÷≈‰ø’º‰
for k=2:1:2*n+1     % Sigma points nolinear propagation
    Y(:,k) = feval(hfx, X(:,k), tpara);
    y = y + Wm(k)*Y(:,k);
end
Pyy = zeros(n); Pxy = zeros(n,m);
for k=1:1:2*n+1
    yerr = Y(:,k)-y;
    Pyy = Pyy + Wc(k)*(yerr*yerr');  % variance
    xerr = X(:,k)-x;
    Pxy = Pxy + Wc(k)*xerr*yerr';  % covariance
end
end