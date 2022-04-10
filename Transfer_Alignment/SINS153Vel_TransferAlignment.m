function res = SINS153Vel_TransferAlignment(imu, mavp, imuerr, davp,avp0)
%% -----------Introduction------------
%速度匹配15维传递对准算法 ，仅估计杆臂
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------davp : 用于设置Kalman P阵 18*1
%-------imuerr : 用于设置Kalman Q阵 
%-------avp0  ； 初始姿态信息  9*1
%-------mavp: 主惯导信息  单位：N*10 弧度 m/s 弧度 m
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%-------res.xkpk: 2*N  估计值与协方差阵
%%  data length and time
N = length(imu(:,end));
%% kalmman参数初始化
kf = []; m  = 6; n =15;
kf.m = m;  kf.n = n;
kf.Pk = 1*diag([davp(1:3); davp(4:6);imuerr.eb; imuerr.db;[2;2;2]])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb;zeros(kf.n-6,1)]).^2;
kf.Gammak = eye(kf.n);
kf.I = eye(kf.n);
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10);
xkpk = zeros(N,2*kf.n);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
factor = 1.0;
timebar(1, N, '姿态速度匹配算法（可以估计安装误差角）.');             % 第一个参数表示间隔数，第二个参数表示总进度
%% Function realize
for i = 2:N        
    wbib = imu(i,1:3)' ;
    fb = imu(i,4:6)' ;
    dt = imu(i,end) - imu(i-1,end);
    t(i) = imu(i,end);    
    % 传递主惯导的速度和位置
    posm = mavp(i,7:9)';
    vnm = mavp(i,4:6)';
%     [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,posm,vnm,dt); 
    eth = EarthParameter(posm,vnm);
    vel = vel + (Cnb*fb + eth.gcc)*dt;
    qua = qupdt(qua,(wbib - Cnb'*eth.wnin)*dt);
    pos = posm;
    % Kalman滤波
    % 一步预测
    Maa=zeros(3,3); 
    Maa(2)=-eth.wnin(3); 
    Maa(3)= eth.wnin(2); 
    Maa(4)= eth.wnin(3); 
    Maa(6)=-eth.wnin(1); 
    Maa(7)=-eth.wnin(2); 
    Maa(8)= eth.wnin(1); 
    Mva=zeros(3,3);
    fn = Cnb*fb;
    Mva(2)= fn(3); 
    Mva(3)=-fn(2); 
    Mva(4)=-fn(3); 
    Mva(6)= fn(1); 
    Mva(7)= fn(2); 
    Mva(8)=-fn(1); 
    Ft = [Maa zeros(3,3) -Cnb       zeros(3,3)  zeros(3,3);        % 
             Mva zeros(3,3) zeros(3,3) Cnb         zeros(3,3);         
             zeros(9,15);];
    Fk = Ft*dt;
    kf.Phikk_1 = kf.I  + Fk;
     kf.Xk = kf.Phikk_1*kf.Xk; 
    kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
    kf.Qk = kf.Qt*dt;
    kf.Pk = kf.Phikk_1*kf.Pk*kf.Phikk_1' + kf.Gammak*kf.Qk*kf.Gammak';        
    % ------------------量测更新-----------------
    Zk = vel-vnm; 
    kf.Xkk_1 = kf.Xk;
    kf.Pkk_1 = kf.Pk;
    web = wbib- Cnb'*eth.wnie;
    kf.Hk = [zeros(3,3) eye(3)     zeros(3,6)  Cnb*askew(web)];
    kf.Rk = diag([0.1; 0.1; 0.1])^2;       
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
    qua  = qdelphi(qua, factor* kf.Xk(1:3));
    kf.Xk(1:3) = (1-factor)*kf.Xk(1:3);  
    att = q2att(qua);
    Cnb = q2mat(qua);
    % 速度
    vel = vel-factor*kf.Xk(4:6);  
    kf.Xk(4:6) = (1-factor)*kf.Xk(4:6);  
    % 记录参数
    avp(i,:) = [att',vel',pos',t(i)];
    xkpk(i,:) = [xk',diag(kf.Pk)'];
    timebar;                                                               % 进度条
end
%% Return Results and Save Workspace
res = varpack(avp,xkpk); 
end

%% End

