function kf_res = kf_phi(imu,gnss,davp,imuerr,avp0)
%% -----------Introduction------------
%15维基本卡尔曼滤波,增加航向观测
%更换时间融合方式
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
L_GNSS = length(gnss(:,end));
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
mbatt = [0;0;0]*pi/180;
Cmb = a2mat(mbatt);%第一次估计出的安装角进行补偿，矩阵是Cmb
% Cmb = (a2mat([-0.0562118148472312,-1.26473527842172e-05,0]))';
ki =1;
gnss(find(gnss(:,7)>7),:) = [];
imugpssyn(imu(:,7), gnss(:,end)); %初始化
gyroStatic = 0.012;accStatic = 0.2;
turnThreshold = 0.004;  %f1数据的设置
% gyroStatic = 7.5e-3;accStatic = 8.5e-3; %仿真数据
% turnThreshold = 3.5e-8; 
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10); 
xkpk = zeros(N ,2*kf.n);   
vm0 = zeros(3,N);
M = zeros(N,6);
D  = zeros(N,6);
TD = zeros(N,1);
mbpitch = zeros(N,1);
mbphi = zeros(N,1);
% meanAtt = zeros(N,2);
zk = zeros(length(gnss),7);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
vm = Cmb*Cnb'*vel;vm0(:,1) = vm;
M(1,1:6) = imu(1,1:6);
D(1,:) = zeros(1,6);
TD(1) = 0;
% j =1;
windowTh = 400;
% meanAtt(1,:) = [mbpitch(1) mbphi(1)];
%% Algorithm develop
timebar(1,N,'SINS/GPS156 Simulation');
for i= 2:N 
        wbib = imu(i,1:3)' ;
        fb = imu(i,4:6)' ;
        dt = imu(i,end) - imu(i-1,end);
        t(i) = imu(i,end);
        %% 
        M(i,:) = (windowTh-1)/windowTh *M(i-1,:) + 1/windowTh *imu(i,1:6);
        D(i,:) = (windowTh-1)/windowTh *D(i-1,:) +  1/windowTh * abs((imu(i,1:6) - M(i,:)));% 静止模式识别
        TD(i) = (windowTh-1)/windowTh *TD(i-1) + 1/windowTh *sum(imu(i,1:3).^2);% 转弯和直线行驶模式识别  
        if D(i,1:3) < gyroStatic & D(i,4:6) < accStatic
            StaticFlag = 0;
        else 
            StaticFlag = 1;% no static
        end
        if TD(i) < turnThreshold
            turnFlag = 1; %no turn
        else
            turnFlag =0;
        end
        %% 惯导更新
        [att,vel,pos,qua,Cnb,eth] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);    
        %-------------kalman预测-------------------      
        kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
        kf.Xk = kf.Phikk_1*kf.Xk; xk = kf.Xk;
        kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
        kf.Qk = kf.Qt*dt;
        kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
        %-------------kalman更新------------------- 		
        [kgps, ~] = imugpssyn(i, i, 'F');
        if kgps>0  %&& kgps <= L_GNSS
            kf.Xkk_1 = kf.Xk;
            kf.Pkk_1 = kf.Pk;
            if sqrt(gnss(kgps,1).^2 + gnss(kgps,2).^2 + gnss(kgps,3).^2) >10 && turnFlag
                phiGPS = atan2(-gnss(kgps,1),gnss(kgps,2));     
                Zk = [att(3) -  phiGPS;%直接做差是否有问题
                      vel-gnss(kgps,1:3)';
                      pos-gnss(kgps,4:6)'];	
                kf.Hk =zeros(7,kf.n);kf.Hk(2:7,4:9) = eye(6);kf.Hk(3,3) = 1;
                kf.Rk = diag([0.2*pi/180;davp(4:6);davp(7:9)]).^2;  
                zk(i,:) = Zk'; 
            else
                Zk = [vel-gnss(kgps,1:3)';
                      pos-gnss(kgps,4:6)'];	
                kf.Hk =zeros(kf.m,kf.n);kf.Hk(1:6,4:9) = eye(6);
                kf.Rk = diag([davp(4:6);davp(7:9)]).^2;
                zk(i,2:7) = Zk';
             end
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
            v = sqrt(vel(1).^2 + vel(2).^2 + vel(3).^2);
%             if (v > 10 ) & turnFlag & diag(kf.Pk(4:6,4:6)) < 0.05
%                  att(3) = atan2(-vel(1),vel(2));% Cnb = a2mat(att);qua = a2qua(att);        
%             end 
        if   v > 10 && turnFlag  &&  StaticFlag
%             att(3) = atan2(-vel(1),vel(2)); %Cnb = a2mat(att);qua = a2qua(att);
            Vb  = Cnb'* vel;
            Vm = sqrt(sum(Vb.^2));
            mbpitch(i) = asin(Vb(3)/Vm);
            mbphi(i) = atan2(-Vb(1),Vb(2));
            Cbm = a2mat([mbpitch(i);0;mbphi(i)]);Cmb =  Cbm';
%             j = j+1;
        else
           mbpitch(i) = mbpitch(i-1) ;
           mbphi(i) = mbphi(i-1);
        end
%         meanAtt(i,:) = meanAtt(i-1,:) +  ([mbpitch(i) mbphi(i)] - meanAtt(i-1,:))/(i+1);
        vm = Cmb*Cnb'*vel;vm0(:,i) = vm;
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
kf_res = varpack(avp,xkpk,zk,vm0,mbpitch,mbphi);
end