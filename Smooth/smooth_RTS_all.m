function  RTS_res = smooth_RTS_all(imu,gnss,davp,imuerr,avp0)
%% -----------Introduction------------
%15维RTS平滑卡尔曼滤波 滑动区间取整个数据时间
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------gnss: 卫星信号 N*8 倒数第二列是定位精度因子  单位： m/s     rad rad m 
%-------davp : 用于设置Kalman P阵 15*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*15 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
L_GNSS = length(gnss(:,end));
%% kalmman参数初始化
kf = []; m  = 6; n =15;
kf.m = m;  kf.n = n;
kf.Pk = 10*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;zeros(kf.n-6,1)]).^2;
kf.Gammak = eye(kf.n);
kf.I = eye(kf.n);
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10);
avp_smooth = zeros(N,10);
dx = zeros(kf.n, N);%量测更新的状态
dx_timeupd = zeros(kf.n, N);%时间更新的状态
dx_smooth = zeros(kf.n, N);%平滑的状态
cov = zeros(kf.n, N);
cov_smooth = zeros(kf.n,N);
P = zeros(kf.n,kf.n,N);
P_timeupd = zeros(kf.n,kf.n,N);
P_smooth = zeros(kf.n,kf.n,N);
F = zeros(kf.n,kf.n,N);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
dx(:,1) = kf.Xk;
dx_timeupd(:,1) =  kf.Xk;

P_timeupd(:,:,1) = kf.Pk;
P(:,:,1)= kf.Pk;
avp(1,:) = [att',vel',pos',t(1)];
avp_smooth(1,:) = avp(1,:);
%% others parameter setting
ki = 1;
%% Algorithm develop
timebar(1,N,'RTS平滑 滑动窗长度为N.');
for i = 2:N 
    wbib = imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i,end) - imu(i-1,end);
    t(i) = imu(i,end);
    %% 惯导更新
    [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
    eth = EarthParameter(pos,vel);%更新当前时刻的曲率半径
    %-------------kalman预测-------------------      
    kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
    kf.Xk = kf.Phikk_1*kf.Xk; xk = kf.Xk;
    kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
    kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    dx_timeupd(:,i) =  xk;%存储时间更新的状态量
    P_timeupd(:,:,i) = kf.Pk;%存储时间更新的状态协方差量
    F(:,:,i) = kf.Phikk_1;
    if (ki<=L_GNSS && gnss(ki,7)> 7)   %根据DPOP判断数据是否可用
        gnss(ki,:) = [];
        L_GNSS =L_GNSS - 1;
    end     
    if ki<=L_GNSS && gnss(ki,end)<=imu(i,end) 
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk;
        Zk = [vel-gnss(ki,1:3)';
              pos-gnss(ki,4:6)'];
        kf.Hk =zeros(kf.m,kf.n);kf.Hk(1:6,4:9) = eye(6);
        kf.Rk = diag([davp(4:6);davp(7:9)]).^2;
        kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
        kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
        kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
        kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
        kf.Pk = (kf.Pk+kf.Pk')/2; 
        [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);
        ki = ki+1;           
    end        
    avp(i,:) = [att',vel',pos',t(i)];
    dx(:,i) = xk;    
    P(:,:,i) = kf.Pk;
    cov(:, i) = diag(kf.Pk);
timebar;  
end
%-----------------后向滤波-------------------------
%% Initial data,利用正向的结果
dx_smooth(:,N)= dx(:,N);
P_smooth(:,:,N)=P(:,:,N);
cov_smooth(:,N) = diag(P_smooth(:,:,N));
avp_smooth(N,:) = avp(N,:);
for i = N-1:-1:1 
    %RTS算法
    Ks_k = P(:,:,i)*F(:,:,i)'*invbc(P_timeupd(:,:,i+1));
    P_smooth(:,:,i) = P(:,:,i) + Ks_k*( P_smooth(:,:,i+1)-P_timeupd(:,:,i+1))*Ks_k';
    dx_smooth(:,i) = dx(:, i) + Ks_k*(dx_smooth(:,i+1)-dx_timeupd(:,i+1)); 
    avp_smooth(i,4:9) = avp(i,4:9) -1*dx_smooth(4:9,i)'; 
    Cn_b = a2mat(avp(i,1:3));
    Cn_n = a2mat(dx_smooth(1:3,i));
    Cnb = Cn_b*(Cn_n)';
    avp_smooth(i,1:3) =  m2att(Cnb);
    P_smooth(:,:,i) = (P_smooth(:,:,i)+P_smooth(:,:,i)')/2;
    cov_smooth(:,i) = diag(P_smooth(:, :, i));
end
RTS_res = varpack( avp_smooth,cov_smooth, P, P_smooth, dx, dx_smooth);
end