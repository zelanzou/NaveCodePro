function res = VPmaster(imu,mavp,davp,imuerr,avp0)
%% -----------Introduction------------
%利用主惯导做速度位置组合，数据量太大，应该将主惯导的频率减小，或者适用部分主惯导信息卡尔曼更新的计算才不会那么大
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------gnss: 卫星信号 N*8 倒数第二列是定位精度因子  单位： m/s     rad rad m 
%-------davp : 用于设置Kalman P阵 15*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
% L_GNSS = length(mavp(:,end));
%% kalmman参数初始化
kf = []; m  = 6; n =15;
kf.m = m;  kf.n = n;
kf.Pk = 10*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;zeros(kf.n-6,1)]).^2;
kf.Gammak = eye(kf.n);
kf.adaptive = 0; 
kf.pconstrain = 0;
kf.I = eye(kf.n);
%% others parameter setting
ki =1;
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10); 
xkpk = zeros(N ,2*kf.n);   
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
%% Algorithm develop
timebar(1,N,'SINS/GPS156 Simulation');
for i= 2:N 
        wbib = imu(i,1:3)' ;
        fb = imu(i,4:6)' ;
        dt = imu(i,end) - imu(i-1,end);
        t(i) = imu(i,end);
        %% 惯导更新
        [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
        eth = EarthParameter(pos,vel);%更新当前时刻的曲率半径
        %-------------kalman预测-------------------      
        kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
        kf.Xk = kf.Phikk_1*kf.Xk; 
        kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
        kf.Qk = kf.Qt*dt;
        kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
        %-------------kalman更新------------------- 		
            kf.Xkk_1 = kf.Xk;
            kf.Pkk_1 = kf.Pk;
            Zk = [vel-mavp(i,4:6)';
                  pos-mavp(i,7:9)'];	
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
            zk(i,:) = Zk';
            ki = ki+1;
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
res = varpack(avp,xkpk);
end