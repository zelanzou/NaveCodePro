function res = UKF_constrait617(imu,davp,imuerr,avp0)
%% -----------Introduction------------
%6*17维车辆约束UKF Kalman滤波  纯惯性下的约束 无GNSS
%观测量：天向速度，高度，车体侧向和天向速度，俯仰姿态角 横滚姿态角
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------davp : 用于设置Kalman P阵 18*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
%% kalmman参数初始化
kf = []; m  = 6; n =17;
kf.m = m;  kf.n = n;
kf.Pk = 10*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db;[5;5]*pi/180])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;]).^2;
kf.Gammak = eye(kf.n);
kf.I = eye(kf.n);
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10);
xkpk = zeros(N ,2*kf.n);
Vu =  zeros(N,1);
h = zeros(N,1);
pitch  = zeros(N,1);
roll = zeros(N,1);
vm0 = zeros(3,N);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
%% others parameter setting
Vu(1) = vel(3);h(1) = pos(3);pitch(1) = att(1);roll(1) = att(2);
mbatt = [0;0;0];
Cmb = a2mat(mbatt.*glv.deg);
[gamma,Wm,Wc] = UKFParameter(kf.n);

ki = timebar(1, N, 'UKF+内部约束模型验证.');
for i = 2:N
    wbib = imu(i,1:3)' ;%当反馈的零偏并不正确时，会导致imu数据不正确，出现波纹状
    fb = imu(i,4:6)' ;
    %% 惯导更新
    [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
    eth = EarthParameter(pos,vel);%更新当前时刻的曲率半径
    Vu(i) = vel(3);h(i) = pos(3);pitch(i) = att(1);roll(i) = att(2);
    Cbn = Cnb';
    vm = Cmb*Cbn*vel;vm0(:,i) = vm;
    M2 = -Cmb*Cbn*askew(vel);M1 = Cmb*Cbn; M3 = Cmb*askew(Cbn*vel);        
    Z = [Vu(i)-Vu(i-1);
          h(i)-h(i-1);
          vm(1);
          vm(3);
          pitch(i)-pitch(i-1);
          roll(i)-roll(i-1);
          ];
    H = [
         zeros(1,5) 1       zeros(1,11) 
         zeros(1,8) 1       zeros(1,8) 
         M2(1,:)    M1(1,:) zeros(1,9) 0       M3(1,3);
         M2(3,:)    M1(3,:) zeros(1,9) M3(3,1) 0;
         1          0       zeros(1,15);
         0          1       zeros(1,15);
        ];
    kf = ukf_filter(kf,[att;vel;pos],imu,gamma,Wm,Wc,nts,1,Z,H);
	xk = kf.Xk;
    %进行反馈
    vel = vel - kf.Xk(4:6);
    pos = pos- kf.Xk(7:9);
    qua = qdelphi(qua, kf.Xk(1:3));Cnb = q2mat(qua);att = m2att(Cnb);
    kf.Xk(1:9)= 0;   
	xkpk(i,:) = [xk',diag(kf.Pk)']; 
    avp(i,:) = [att',vel',pos',t(i)];
    timebar;
end 		
%% Save Workspace
res =  varpack(avp,xkpk); 
end
