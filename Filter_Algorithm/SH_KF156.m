function kf_res = SH_KF156(imu,gnss,davp,imuerr,avp0)
%% -----------Introduction------------
%Sage-Husa自适应卡尔曼滤波
%需要注意的问题：1.同步时间必须处理好，可以用imugpsavpr数据测试 2.F1数据误差增大不正常
%自适应量测与方差受限
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------gnss: 卫星信号 N*7 包含定位精度因子  单位： m/s     rad rad m 
%-------davp : 用于设置Kalman P阵 15*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
L_GNSS = length(gnss(:,end));
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
%% kalmman参数初始化
kf = []; m  = 6; n =15;
kf.m = m;  kf.n = n;
kf.Pk = 1*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;zeros(kf.n-6,1)]).^2;
kf.Rk = diag([davp(4:6);davp(7:9)]).^2;
kf.Gammak = eye(kf.n); kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
kf.I = eye(kf.n);
%% others parameter setting
kf.b = 0.92;kf.beta = 1;
kf.Rmin = diag(diag(kf.Rk)/200);
kf.Rmax = 200*kf.Rk;
ki =1;
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);

t(1) = imu(1,end);
avp = zeros(N,10); avp(1,:) = [att',vel',pos',t(1)];
xkpk = zeros(N ,2*kf.n);   xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
%% Algorithm develop
timebar(1,N,'15维量测自适应+方差受限');
for i= 2:N 
        wbib = (imu(i,1:3)'+imu(i-1,1:3)')/2 ;
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
        if ki<=L_GNSS && gnss(ki,end)<=imu(i,end)  
            kf.Xkk_1 = kf.Xk;
            kf.Pkk_1 = kf.Pk;
            Zk = [vel-gnss(ki,1:3)';
                     pos-gnss(ki,4:6)'];
            kf.Hk =zeros(kf.m,kf.n);kf.Hk(1:6,4:9) = eye(6);
            kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
            kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
			kf.Py0 =  kf.Hk*kf.PXZkk_1;
            for k=1:kf.m
                ry = kf.rk(k)^2 -  kf.Py0(k,k); %序贯量测
                if ry<kf.Rmin(k,k)
                    kf.Rk(k,k)  = (1-kf.beta)*kf.Rk(k,k) +kf.beta*kf.Rmin(k,k); 
                elseif ry>kf.Rmax(k,k)
                    kf.Rk(k,k) =  kf.Rmax(k,k);
                else
                    kf.Rk(k,k) = (1-kf.beta)*kf.Rk(k,k) + kf.beta*ry;
                end
            end
            kf.beta = kf.beta/(kf.beta+kf.b);	
                
            kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
            kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
            kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
            kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
            kf.Pk = (kf.Pk+kf.Pk')/2;
 			
            [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua); 
            zk(ki,:) = [Zk',sqrt(diag(kf.Rk))',ki];
            ki = ki+1;
        end
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
kf_res = varpack(avp,xkpk,zk);
end