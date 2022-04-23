function res = AttAddKmtcCons(imu, mavp, imuerr, davp,avp0,index)
%% -----------Introduction------------
%姿态匹配+ 车体运动学速度匹配 20维传递对准算法 估计主子安装角，车体安装角，不能估计杆臂
%可以分析gnss缺失的情况
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------davp : 用于设置Kalman P阵 18*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%-------mavp: 主惯导信息  单位：N*10 弧度 m/s 弧度 m
%-------index: gnss缺失区间
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
%% kalmman参数初始化
kf = []; m  = 5; n =20;
kf.m = m;  kf.n = n;
kf.Pk = 10*diag([davp(1:3); davp(4:6); davp(7:9);imuerr.eb; imuerr.db;[10;10;10]*pi/180;[10;10]*pi/180])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;zeros(kf.n-6,1)]).^2;
kf.Gammak = eye(kf.n);
kf.I = eye(kf.n);
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10);
xkpk = zeros(N,2*kf.n);
vm0 = zeros(3,N);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
%% others parameter setting
ftrsn = 1;% 传递对准标志位
fpk = 0;
factor = 1.0;
mbatt = [0;0;0]*pi/180;
Cmb = a2mat(mbatt);
k=1;
%% Algorithm develop
timebar(1,N, '姿态匹配与运动学约束多源融合算法.');
for i = 2:N    
    wbib = imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i,end) - imu(i-1,end);
    t(i) = imu(i,end);  
    if k<length(index) && i>index(k)*200 && i<index(k+1)*200
        ftrsn = 0;
        fpk = 1;
        k=k+1;
    else
        ftrsn = 1;
        % 传递主惯导的速度和位置
%        pos = mavp(i,7:9)';
%         vel = mavp(i,4:6)';
        vnm = mavp(i,4:6)';
    end
    [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt); 
    eth = EarthParameter(pos,vel); 
    Cbn = Cnb';
    vm = Cmb *Cbn*vel;vm0(:,i) = vm;
    %-------------kalman预测-------------------      
    kf.Phikk_1 = kf.I + KF_Phi(eth,Cnb,fb,kf.n)*dt;%离散化二阶泰勒展开
    kf.Xk = kf.Phikk_1*kf.Xk; 
    kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
    kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';   
    % 量测更新
    kf.Xkk_1 = kf.Xk;
    kf.Pkk_1 = kf.Pk;
    if ftrsn == 1        
        % 第二次进入，要重新修改pk
        if fpk == 1
             kf.Pkk_1 = 10* kf.Pkk_1;
             fpk =0;
        end        
        % 计算姿态误差
        mqnb = a2qua(mavp(i,1:3));
        sqbn = qconj(qua);
        Cnbm = q2mat(mqnb);
        Cbsn = q2mat(sqbn);
        z = Cnbm*Cbsn; 

        M2 = -Cmb*Cbn*askew(vel);M1 = Cmb*Cbn;M3 = askew(vm);     
        kf.Hk = [eye(3)  zeros(3,3) zeros(3,9) -Cnb     zeros(3,2);                                             % 
                 M2(1,:) M1(1,:)    zeros(1,12) 0       M3(1,3);
                 M2(3,:) M1(3,:)    zeros(1,12) M3(3,1) 0;
                 ];
        kf.Rk = diag([[100; 100; 100;].*pi/180/60;[0.5;0.5]])^2;                      % 姿态速度匹配  
        Zk = [0.5*(z(3,2)-z(2,3)); 
              0.5*(z(1,3)-z(3,1)); 
              0.5*(z(2,1)-z(1,2));                                             % 姿态
              vm(1);
              vm(3);];     
    else         
        M2 = -Cmb*Cbn*askew(vel);M1 = Cmb*Cbn;M3 = askew(vm); 
        kf.Hk = [M2(1,:)  M1(1,:) zeros(1,12) 0       M3(1,3);
                   M2(3,:)  M1(3,:) zeros(1,12) M3(3,1) 0;]; 
        kf.Rk = diag([0.8;0.8])^2;     
        Zk = [vm(1);vm(3)];                                              % 速度量测
    end
    kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
    kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
    kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
    kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
    kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
    kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
    kf.Pk = (kf.Pk+kf.Pk')/2;
    % 反馈 
    xk = kf.Xk;
    % 姿态
    qua  = qdelphi(qua, factor*kf.Xk(1:3));
    kf.Xk(1:3) = (1-factor)*kf.Xk(1:3);  
    att = q2att(qua);
    Cnb = q2mat(qua);
    % 速度
    vel = vel-factor*kf.Xk(4:6);  
    kf.Xk(4:6) = (1-factor)*kf.Xk(4:6);  
    %位置
    if ftrsn == 1  
        pos = mavp(i,7:9)';
        kf.Xk(7:9)  = zeros(3,1);
    else
        pos =  pos - factor*kf.Xk(7:9);  
        kf.Xk(7:9) = (1-factor)*kf.Xk(7:9); 
    end
    % 记录参数
    avp(i,:) = [att',vel',pos',t(i)];
    xkpk(i,:) = [xk',diag(kf.Pk)'];
    timebar;   
end 
%% 返回结果
res = varpack(xkpk, avp, vm0); 
% save 'RES.mat'
end

%% 函数结束