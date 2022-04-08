function [att,p,v,qua,Cnb,kf,xk]= feed_back_correct(kf,avp,qua)
%·´À¡Ð£Õý
ratio = 0.8;
v = avp(4:6) - ratio *kf.Xk(4:6);
p = avp(7:9) - ratio *kf.Xk(7:9);
qua = qdelphi(qua, ratio*kf.Xk(1:3));Cnb = q2mat(qua);
% dq = rv2q(ratio*kf.Xk(1:3));Cnb_n1 = q2mat(dq);Cnbk=Cnb_n1*Cnbk;
% qua = m2qua(Cnbk);
% qua =  qmul(rv2q(ratio *kf.Xk(1:3)), qua);
% Cnbk = q2mat(qua);
att = m2att(Cnb);
% kf.Xk(1:9) = 0;
xk = kf.Xk;
kf.Xk(1:9)= (1-ratio)*kf.Xk(1:9);
end