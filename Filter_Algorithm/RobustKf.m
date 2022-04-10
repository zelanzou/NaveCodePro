function kf_res = RobustKf(imu,gnss,davp,imuerr,avp0)
%ref:[1].A Novel Adaptively-Robust Strategy Based on the Mahalanobis Distance for GPS/INS Integrated Navigation Systems
%    [2].Robust Kalman filtering based on Mahalanobis distance as outlier judging criterion
%% -----------Introduction------------
%15维基本卡尔曼滤波
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
ki =1;
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);
avp = zeros(N,10); 
xkpk = zeros(N ,2*kf.n);   
lamda = zeros(6,1);
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
        [att,vel,pos,qua,Cnb,eth] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
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
%             zk(ki,:) = kf.rk';                  
            kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
            kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%;理论残差协方差
%             thecov = kf.PZZkk_1;
%             estcov = kf.rk*kf.rk';   %实际协方差
%             alpha(ki) = kf.rk'*kf.rk/trace(kf.PZZkk_1); %自适应因子1
            
            Tk = kf.rk'*inv(kf.PZZkk_1)*kf.rk;
%             A(ki) = Tk;
            if Tk >18.55 %显著性水平0.005
                lamdak = Tk/18.55;
            else
                lamdak = 1;
            end
            
            strk(:,1) = abs(kf.rk)./sqrt(diag(kf.PZZkk_1)); %标准化残差
            %使用IGIII降权函数，根据标准化残差构建权因子
            for g= 1:6
                if strk(g,1) <2
                    lamda(g) = 1;
                elseif strk(g,1) >2 && strk(g,1) <4.5
                     lamda(g) = 2/strk(g,1)*((4.5-strk(g,1))/(4.5-2))^2;
                else
                    lamda(g) = 1e-8;
                end
            end
            
%             l(ki,:) = 1./lamda;
%             la(ki,:) = lamdak;
            kf.Rk = 1./lamda.*kf.Rk; %等价权矩阵
%             R(ki,:) = diag(kf.Rk)';
            kf.PZZkk_1 = (lamdak)*kf.Hk*kf.PXZkk_1 + kf.Rk;%;重新计算
            
            kf.Kk = kf.PXZkk_1*(1/lamdak)*invbc(kf.PZZkk_1);
            kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
            kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
            kf.Pk = (kf.Pk+kf.Pk')/2;
            
            [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);             
            ki = ki+1;
        end
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
kf_res = varpack(avp,xkpk);
end