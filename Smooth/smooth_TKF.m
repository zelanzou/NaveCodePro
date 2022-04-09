function TKF_res = smooth_TKF(imu,gnss,davp,imuerr,avp0)
%% -----------Introduction------------
%15维双向平滑卡尔曼滤波
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
%% others parameter setting
ki =1;
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,kf.n+1);  %便于融合，将零偏估计同时存储
dx_mea_upd = zeros(N,kf.n);
cov_p_mea = zeros(N,kf.n); 
% xkpk = zeros(N ,2*kf.n);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
dx_mea_upd(1,:) = kf.Xk;
cov_p_mea(1,:) = diag(kf.Pk);
gb_est = kf.Xk(10:12);
ab_est = kf.Xk(13:15);
avp(1,:) = [att',vel',pos',gb_est',ab_est',t(1)];
%% Algorithm develop
timebar(1,N,'forward processing.');
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
    kf.Xk = kf.Phikk_1*kf.Xk; xk = kf.Xk;
    kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
    kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    %-------------kalman更新------------------- 		
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
    dx_mea_upd(i,:) = xk;
    cov_p_mea(i,:) = diag(kf.Pk);
    gb_est = kf.Xk(10:12);
    ab_est = kf.Xk(13:15);
    avp(i,:) = [att',vel',pos',gb_est',ab_est',t(i)];
    timebar;
end
%---------------------逆向滤波------------------------
%% others parameter setting
vel = -vel;
ikf = kf;
idx = [4:6,10:12];
ikf.Pk = 10*ikf.Pk; %反向P阵设大
%% Memory allocation
iavp = zeros(N,ikf.n+1);% 
idx_mea_upd = zeros(N,ikf.n);
icov_p_mea = zeros(N,ikf.n);
%% Initial data
ikf.Xk(idx) = -ikf.Xk(idx);ikf.Xk(1:9) = 1.5*ikf.Xk(1:9);
iavp(N,:) = avp(end,:);iavp(N,idx) = -iavp(N,idx);
idx_mea_upd(N,:) = ikf.Xk;
icov_p_mea(N,:) = diag(ikf.Pk);
ki = L_GNSS;
%% Algorithm develop
timebar(1,N,'backward processing.');
for i = N-1:-1:1
    wbib = -imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i+1,end) - imu(i,end);
    t(i) = imu(i,end);
%---------------纯惯性解算    
    ieth = reverseEarthParameter(pos,vel);
    wbnb = wbib - Cnb'*ieth.wnin;%b系下
    qua = qupdt(qua,wbnb*dt);%wbnb*dt是相当于等效旋转增量rv
    Cnb = q2mat(qua);
    att = m2att(Cnb);

    fn = Cnb*fb;
    an = fn + ieth.gcc;
    vel = vel+ an*dt;

    Mpv = [0 1/ieth.RMh 0;1/ieth.clRNh 0 0;0 0 1];
    pos = pos + Mpv*vel*dt;
%     Cnb = Cnbk;   
    
    ieth = reverseEarthParameter(pos,vel);%更新当前时刻的曲率半径
    ikf.Phikk_1 = ikf.I + KF_Phi(ieth,Cnb,fb,15)*dt;%离散化二阶泰勒展开
    ikf.Xk = ikf.Phikk_1*ikf.Xk;xk = ikf.Xk;
    ikf.Gammak(1:3,1:3) = -Cnb; ikf.Gammak(4:6,4:6) = Cnb;
    ikf.Qk = ikf.Qt*dt;
    ikf.Pk = ikf.Phikk_1*ikf.Pk*ikf.Phikk_1' + ikf.Gammak*ikf.Qk*ikf.Gammak';
     if ( (ki<=L_GNSS)&& (imu(i,end)-gnss(ki,end)<=1e-6))   %防止matlab数据格式不一致
        Zk = [vel-(-gnss(ki,1:3))';
              pos-gnss(ki,4:6)'];	
        ikf.Xkk_1 = ikf.Xk;
        ikf.Pkk_1 = ikf.Pk;
        ikf.PXZkk_1 = ikf.Pkk_1*ikf.Hk';%状态一步预测与量测一步预测的协均方误差	
        ikf.rk = Zk - ikf.Hk*ikf.Xkk_1;%残差        
        ikf.PZZkk_1 = ikf.Hk*ikf.PXZkk_1 + ikf.Rk;%量测一步预测均方误差阵
        ikf.Kk = ikf.PXZkk_1*invbc(ikf.PZZkk_1);
        ikf.Xk = ikf.Xkk_1 + ikf.Kk*ikf.rk;
        ikf.Pk = ikf.Pkk_1 - ikf.Kk*ikf.PZZkk_1*ikf.Kk';
        ikf.Pk = (ikf.Pk+ikf.Pk')/2;
        [att,pos,vel,qua,Cnb,ikf,xk]= feed_back_correct(ikf,[att;vel;pos],qua);
        ki = ki-1;
        if ki == 0
            ki=1;
        end
     end
    idx_mea_upd(i,:) = xk;
    icov_p_mea(i,:) = diag(ikf.Pk);
    gb_est = ikf.Xk(10:12);
    ab_est = ikf.Xk(13:15);    
    iavp(i,:) = [att',vel',pos',gb_est',ab_est',t(i)];
    timebar;
end
iavp(:,idx) = -iavp(:,idx);
ps = cov_p_mea + icov_p_mea;  % 严教材 P149  公式6.3.16  注意：融合的是avp结果，而不是状态参数的误差
xs = icov_p_mea./(cov_p_mea + icov_p_mea).*avp(:,1:end-1)+...
cov_p_mea./(cov_p_mea + icov_p_mea).*iavp(:,1:end-1);

TKF_res = varpack(avp,iavp,xs,ps,cov_p_mea,dx_mea_upd,icov_p_mea,idx_mea_upd);
end