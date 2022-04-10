function kf_res = P_AdaptiveKf_NHC_PR(imu,gnss,davp,imuerr,avp0,Rmin,Rmax)
%% -----------Introduction------------
%17维位置组合卡尔曼滤波 ，考虑GNSS缺失情况下使用运动学约束 , 模式识别（PR）
%解析计算车体安装角
%车体安装角只在约束阶段进行估计
%车体约束以SINS频率更新
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------gnss: 卫星信号 N*8 倒数第二列是定位精度因子  单位： m/s     rad rad m 
%-------davp : 用于设置Kalman P阵 17*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
%%  data length and time
N = length(imu(:,end));
L_GNSS = length(gnss(:,end));
%% kalmman参数初始化
kf = []; m  = 3; n =15;
kf.m = m;  kf.n = n;
kf.Pk = 10*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db;])^2;%[5;5]*pi/180
kf.Qt = diag([imuerr.web; imuerr.wdb; zeros(kf.n-6,1)]).^2;
kf.Gammak = eye(kf.n); 
kf.I = eye(kf.n);

%% constrait parameter setting
upTolow = 1;%高精度到低精度切换滤波器的次数 调小P阵
lowToup = 0;%低精度到高精度切换滤波器的次数 调大P阵
ki =1;
j=1;
mbatt = [0;0;0]*pi/180;
Cmb = a2mat(mbatt);
kf.adaptive = 1; 
kf.b = 0.5;kf.beta = 1;



% kf.pconstrain = 1;
% kf.Pkmax = (diag(kf.Pk)+1)*1e10; 
% kf.Pkmin = Pmin;
% gyroStatic = 0.012;accStatic = 0.2;
turnThreshold = 0.004;
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10); 
xkpk = zeros(N ,2*kf.n);   
vm0 = zeros(3,N);
mbpitch = zeros(N,1);
mbphi = zeros(N,1);
M = zeros(N,6);
D  = zeros(N,6);
TD = zeros(N,1);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
M(1,1:6) = imu(1,1:6);
D(1,:) = zeros(1,6);
TD(1) = 0;
% j =1;
windowTh = 400;
%% Algorithm develop
timebar(1,N,'位置组合+运动学约束');
for i= 2:N 
    wbib = imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i,end) - imu(i-1,end);
    t(i) = imu(i,end);
    %% 
    M(i,:) = (windowTh-1)/windowTh *M(i-1,:) + 1/windowTh *imu(i,1:6);
    D(i,:) = (windowTh-1)/windowTh *D(i-1,:) +  1/windowTh * abs((imu(i,1:6) - M(i,:)));% 静止模式识别
    TD(i) = (windowTh-1)/windowTh *TD(i-1) + 1/windowTh *sum(imu(i,1:3).^2);% 转弯和直线行驶模式识别  
%     if D(i,1:3) < gyroStatic & D(i,4:6) < accStatic
%         StaticFlag = 0;
%     else 
%         StaticFlag = 1;% no static
%     end
    if TD(i) < turnThreshold
        turnFlag = 1; %no turn
    else
        turnFlag =0;
    end
    %% 惯导更新
    [att,vel,pos,qua,Cnb,eth] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
     Cbn = Cnb';
     vm = Cmb *Cbn*vel;vm0(:,i) = vm;
    %% kalman滤波  
	kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
	kf.Xk = kf.Phikk_1*kf.Xk; xk = kf.Xk;
	kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
	kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';
    %% 
    if (ki<=L_GNSS && gnss(ki,7)> 7)   %根据DPOP判断数据是否可用
        gnss(ki,:) = [];
        L_GNSS =L_GNSS - 1;
    end  
    if ki<=L_GNSS && gnss(ki,end)<=imu(i,end)
        if(lowToup == 1)
           upTolow = 1;
           lowToup = 0;
           kf.Pk  = 100* kf.Pk;
        end        
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk;
        Zk = [pos-gnss(ki,4:6)'];            
        kf.Hk =zeros(3,kf.n);kf.Hk(1:3,7:9) = eye(3);
        kf.Rk = diag([davp(7:9)]).^2;
        kf.Rmin = diag(Rmin);
        kf.Rmax = diag(Rmax);
        kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
        kf.Py0 =  kf.Hk*kf.PXZkk_1;
        if kf.adaptive==1  % for adaptive KF, make sure Rk is diag 24/04/2015 ？?
            for k=1:kf.m
                ry = kf.rk(k)^2 -  kf.Py0(k,k); %序贯量测+方差受限
                if ry<kf.Rmin(k,k), ry = kf.Rmin(k,k); end
                if ry>kf.Rmax(k,k),     kf.Rk(k,k) = kf.Rmax(k,k);
                else                	kf.Rk(k,k) = (1-kf.beta)*kf.Rk(k,k) + kf.beta*ry;
                end
            end
            kf.beta = kf.beta/(kf.beta+kf.b);
        end 
        kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
        kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
        kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
        kf.Pk = (kf.Pk+kf.Pk')/2;         
        [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);             
        zpk(ki,:) = [Zk',diag(kf.Rk)'];
        ki = ki+1;
        meaflag =1;
    end
    if  ki >1 &&  ki < length(gnss) && gnss(ki,end)-gnss(ki-1,end) > 1    %已掉帧
        if(upTolow == 1)
            lowToup = 1;
            upTolow = 0;
            kf.Pk  = 1*kf.Pk;
        end        
        M2 = -Cmb*Cbn*askew(vel);M1 = Cmb*Cbn;M3 = -askew(vm);
        Zk = [vm(1);
              vm(3);
                ]; 
        kf.Hk = [M2(1,:)    M1(1,:) zeros(1,9) ;%0       M3(1,3);
                 M2(3,:)    M1(3,:) zeros(1,9) ;%M3(3,1) 0;
                 ];	
        kf.Rk = diag([0.8;0.8]).^2;
        kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
        
        kf.Xkk_1 = kf.Xk;
        kf.Pkk_1 = kf.Pk;
        kf.Rmin = diag([0.05;0.05]).^2;
        kf.Rmax = diag([8;8]).^2;
        kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差        
        kf.Py0 =  kf.Hk*kf.PXZkk_1;
        if kf.adaptive==1  % for adaptive KF, make sure Rk is diag 24/04/2015 ？?
            for k=1:2
                ry = kf.rk(k)^2 -  kf.Py0(k,k); %序贯量测+方差受限
                if ry<kf.Rmin(k,k), ry = kf.Rmin(k,k); end
                if ry>kf.Rmax(k,k),     kf.Rk(k,k) = kf.Rmax(k,k);
                else                	kf.Rk(k,k) = (1-kf.beta)*kf.Rk(k,k) + kf.beta*ry;
                end
            end
            kf.beta = kf.beta/(kf.beta+kf.b);
        end 
        kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
        kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
        kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
        kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
        kf.Pk = (kf.Pk+kf.Pk')/2;          
        [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);
        zk(j,:) = [Zk',diag(kf.Rk)'];
        j =j+1;
        meaflag =0;
    end
    v = sqrt(vel(1).^2 + vel(2).^2 + vel(3).^2);
    if (v > 10 ) && turnFlag &&  meaflag
%         att(3) = atan2(-vel(1),vel(2)); Cnb = a2mat(att);qua = a2qua(att);
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
    avp(i,:) = [att',vel',pos',t(i)];
    xkpk(i,:) = [xk',diag(kf.Pk)'];
    timebar;
end
kf_res = varpack(avp,xkpk,vm0,mbpitch, mbphi);
end