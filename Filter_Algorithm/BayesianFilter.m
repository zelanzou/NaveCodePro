function kf_res = BayesianFilter(imu,gnss,davp,imuerr,avp0)
%% -----------Introduction------------
%15维鲁棒卡尔曼滤波
%可以解决连续野值问题
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------gnss: 卫星信号 N*8 ,倒数第二列是定位精度因子  单位： m/s     rad rad m 
%-------davp : 用于设置Kalman P阵 15*1
%-------imuerr : 用于设置Kalman Q阵 ,imuerr.eb,imuerr.db,imuerr.web,imuerr.wdb
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
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
kf.Rk1 =  diag([davp(4:6);davp(7:9)]).^2;%将单位一致统一为米
kf.Rk2 =  25*diag([davp(4:6);davp(7:9)]).^2;
kf.Hk =zeros(kf.m,kf.n);kf.Hk(1:6,4:9) = eye(6);
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
%% others parameter setting
ki =1;
a1= 0.9;a2 = 0.1;% 在滤波出现异常量测的先验概率
%% Algorithm develop
timebar(1,N,'15维RKF组合导航');
for i= 2:N 
    wbib = imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i,end) - imu(i-1,end);
    t(i) = imu(i,end);
    %% 惯导更新
    [att,vel,pos,qua,Cnb,eth] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
    %-------------kalman预测-------------------      
    kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
    kf.Xk = kf.Phikk_1*kf.Xk; xk = kf.Xk;
    kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
    kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    %-------------kalman更新------------------- 	
    if ki<=L_GNSS && gnss(ki,end)<=imu(i,end)  
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk;
        Zk = [vel-gnss(ki,1:3)';
              pos-gnss(ki,4:6)'];
        Zkk_1 = kf.Hk*kf.Xkk_1; 
        kf.rk = Zk - Zkk_1;%残差 
        E1 = kf.Hk*kf.Pkk_1*kf.Hk'+kf.Rk1;iE1 = inv(E1);
        E2 = kf.Hk*kf.Pkk_1*kf.Hk'+kf.Rk2;iE2 = inv(E2);
        A1 = 1/(1+a2/a1*sqrt(det(E1)/det(E2))*exp(0.5*kf.rk'*(iE1-iE2)*kf.rk));
        A2 = 1-A1;
        kf.Kk = kf.Pkk_1*kf.Hk'*(A1*iE1+A2*iE2);
        kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
        Fk = A1*A2*(iE1-iE2)*kf.rk*kf.rk'*(iE1-iE2)';
        Bk = A1*iE1+A2*iE2-Fk;
        kf.Pk = (eye(15)-kf.Pkk_1*kf.Hk'*Bk*kf.Hk)*kf.Pkk_1;
        [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);             
%         AA(ki) = A1;BB(ki) = A2;%记录概率
        ki = ki+1;
    end
      
    avp(i,:) = [att',vel',pos',t(i)];
    xkpk(i,:) = [xk',diag(kf.Pk)'];
    timebar;
end
kf_res = varpack(avp,xkpk);
end