function  X = state_function(X,imu,avp,nts)
eth =  EarthParameter(avp(7:9),avp(4:6));
Cnb = a2mat(avp(1:3));
fb = imu(4:6)';
flag = 15;
F = KF_Phi(eth,Cnb,fb,flag);
X = X + F*nts*X;
%% 对比下面模型，该模型存在错误
%     Re 		= 6378245;   %地球长半径
%     e 		= 1/298.257;  %地球扁率
%     wie 	= 7.292e-5;  %地球自转角速度
%     afa = X(1:3);gyrobias = X(10:12);accbias = X(13:15);
%     Ve 		= avp(4);
%     Vn 		= avp(5);
%     Vu 		= avp(6);
%     L 		= avp(7);
%     h 		= avp(9);
% 	s1 = sin(afa(1));   s2 = sin(afa(2)); s3 = sin(afa(3));
% 	c1 = cos(afa(1));   c2 = cos(afa(2)); c3 = cos(afa(3));
%     t1 = s1/c1;
%     C_w = [	c2     0  s2;
%             t1*s2  1  -t1*c2;
%             -s2/c1 0  c2/c1];
%     Cn_n = (a2mat(afa))';
%     Cn_b = a2mat(avp(1:3));
%     Fn = Cn_b*imu(4:6)';
%     Rm 		= Re*(1-2*e+3*e*sin(L)^2);
%     Rn 		= Re*(1-e*sin(L)^2);
%     Mav = [0,               -1/(Rm + h), 0;
%            1/(Rn + h),      0,           0;
%            tan(L)/(Rn + h), 0,           0];
%     Map = [0,                                   0, Vn/(Rm + h)^2;
%            -wie*sin(L),                         0, -Ve/(Rn + h)^2;
%            wie*cos(L) + Ve/(cos(L)^2*(Rn + h)), 0, -(Ve*tan(L))/(Rn + h)^2];
%     wnin = [0;wie*cos(L);wie*sin(L)]+[-Vn/(Rm+h);Ve/(Rn+h);Ve*tan(L)/(Rn+h)];
%     sigma_wnin = Mav*X(4:6) + Map *X(7:9);
%     %姿态模型
%     X(1:3) = X(1:3) + nts*C_w*((eye(3) - Cn_n)*wnin + Cn_n*sigma_wnin - Cn_b*gyrobias);
%     
%     Mvv = [(Vn*tan(L))/(Rn + h) - Vu/(Rn + h),       2*wie*sin(L) + (Ve*tan(L))/(Rn + h), - Ve/(Rn + h) - 2*wie*cos(L);
%             - 2*wie*sin(L) - (2*Ve*tan(L))/(Rn + h), -Vu/(Rm + h),                        -Vn/(Rm + h);
%             (2*Ve)/(Rn + h) + 2*wie*cos(L),          (2*Vn)/(Rm + h),                     0];
%     Mvp = [Vn*(2*wie*cos(L) + Ve/(cos(L)^2*(Rn + h))) + 2*Vu*wie*sin(L), 0, (Ve*Vu)/(Rn + h)^2 - (Ve*Vn*tan(L))/(Rn + h)^2;
%           -Ve*(2*wie*cos(L) + Ve/(cos(L)^2*(Rn + h))),                   0, (Ve^2*tan(L))/(Rn + h)^2 + (Vn*Vu)/(Rm + h)^2;
%           -2*Ve*wie*sin(L),                                              0, - Ve^2/(Rn + h)^2 - Vn^2/(Rm + h)^2];
%    %速度模型 
%    X(4:6) = X(4:6) + nts*((eye(3) - Cn_n)*Fn + Mvv*X(4:6) +  Mvp*X(7:9)+ Cn_n'*Cn_b*(accbias));
%    
%    Mpos =[0,0,0,0,                  1/(Rm + h),0,0,                             0,-Vn/(Rm + h)^2;
%           0,0,0,1/(cos(L)*(Rn + h)),0,         0,(Ve*tan(L))/(cos(L)*(Rn + h)), 0,-Ve/(cos(L)*(Rn + h)^2);
%           0,0,0,0,                  0,         1,0,                             0,0];
%   %位置误差方程与线性相同
%   X(7:9) = X(7:9) + nts*(Mpos*X(1:9));   
%% 
% 
%     Re 		= 6378245;   %地球长半径
%     e 		= 1/298.257;  %地球扁率
%     wie 	= 7.292e-5;  %地球自转角速度
%     PhiE = X(1);
%     PhiN = X(2);
%     PhiU = X(3);
%     Ve 		= avp(4);
%     Vn 		= avp(5);
%     Vu 		= avp(6);
%     L 		= avp(7);
%     h 		= avp(9);
%     Cn_b = a2mat(avp(1:3));
%     Fn = Cn_b*imu(4:6)';
%     fe 		= Fn(1);
%     fn 		= Fn(2);
%     fu 		= Fn(3);
%     Rm 		= Re*(1-2*e+3*e*sin(L)^2);
%     Rn 		= Re*(1-e*sin(L)^2);
%     %连续系统状态转换阵 F 的时间更新
%     F    = zeros(15,15);
%     %模型有区别Maa，Mva如何推导出来的？
%     FN   = [                                                  0,                              wie*sin(L) + (Ve*tan(L))/(Rn + h),                                                 cos(PhiU)*(Ve/(Rn + h) + wie*cos(L)) - (Vn*sin(PhiU))/(Rm + h),                                       0,                         -1/(Rm + h),                            0,                                                            0, 0,                                  Vn/(Rm + h)^2;
%                             - wie*sin(L) - (Ve*tan(L))/(Rn + h),                                                              0,                                                 sin(PhiU)*(Ve/(Rn + h) + wie*cos(L)) + (Vn*cos(PhiU))/(Rm + h),                              1/(Rn + h),                                   0,                            0,                                                  -wie*sin(L), 0,                                 -Ve/(Rn + h)^2;
%  cos(PhiU)*(Ve/(Rn + h) + wie*cos(L)) - (Vn*sin(PhiU))/(Rm + h), sin(PhiU)*(Ve/(Rn + h) + wie*cos(L)) + (Vn*cos(PhiU))/(Rm + h), (Ve/(Rn + h) + wie*cos(L))*(PhiN*cos(PhiU) - PhiE*sin(PhiU)) - (Vn*(PhiE*cos(PhiU) + PhiN*sin(PhiU)))/(Rm + h),                         tan(L)/(Rn + h),                                   0,                            0,                          wie*cos(L) + Ve/(cos(L)^2*(Rn + h)), 0,                        -(Ve*tan(L))/(Rn + h)^2;
%                                                               0,                                                            -fu,                                                                                  - fn*cos(PhiU) - fe*sin(PhiU),      (Vn*tan(L))/(Rn + h) - Vu/(Rn + h), 2*wie*sin(L) + (Ve*tan(L))/(Rn + h), - Ve/(Rn + h) - 2*wie*cos(L), Vn*(2*wie*cos(L) + Ve/(cos(L)^2*(Rn + h))) + 2*Vu*wie*sin(L), 0, (Ve*Vu)/(Rn + h)^2 - (Ve*Vn*tan(L))/(Rn + h)^2;
%                                                              fu,                                                              0,                                                                                    fe*cos(PhiU) - fn*sin(PhiU), - 2*wie*sin(L) - (2*Ve*tan(L))/(Rn + h),                        -Vu/(Rm + h),                 -Vn/(Rm + h),                  -Ve*(2*wie*cos(L) + Ve/(cos(L)^2*(Rn + h))), 0,  (Ve^2*tan(L))/(Rn + h)^2 + (Vn*Vu)/(Rm + h)^2;
%                                   - fn*cos(PhiU) - fe*sin(PhiU),                                    fe*cos(PhiU) - fn*sin(PhiU),                                  - fe*(PhiE*cos(PhiU) + PhiN*sin(PhiU)) - fn*(PhiN*cos(PhiU) - PhiE*sin(PhiU)),          (2*Ve)/(Rn + h) + 2*wie*cos(L),                     (2*Vn)/(Rm + h),                            0,                                             -2*Ve*wie*sin(L), 0,            - Ve^2/(Rn + h)^2 - Vn^2/(Rm + h)^2;
%                                                               0,                                                              0,                                                                                                              0,                                       0,                          1/(Rm + h),                            0,                                                            0, 0,                                 -Vn/(Rm + h)^2;
%                                                               0,                                                              0,                                                                                                              0,                     1/(cos(L)*(Rn + h)),                                   0,                            0,                                (Ve*tan(L))/(cos(L)*(Rn + h)), 0,                        -Ve/(cos(L)*(Rn + h)^2);
%                                                               0,                                                              0,                                                                                                              0,                                       0,                                   0,                            1,                                                            0, 0,                                              0 ] ;
% 
%     F(1:9, 1:9)  = FN;                                                       
%     F(1:3,10:12) = Cn_b;%这里没有负号 ？
%     F(4:6,13:15) = Cn_b; 
%     kf.Gammak(1:3,1:3) = -Cn_b; kf.Gammak(4:6,4:6) = Cn_b;
%     A 	 = eye(15,15)+F*nts;
%     X = A*X;
end