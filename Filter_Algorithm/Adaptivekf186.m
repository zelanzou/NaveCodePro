function kf_res = Adaptivekf186(imu,gnss,davp,imuerr,avp0,Rmin,Pmin)
%% -----------Introduction------------
%18维自适应卡尔曼滤波
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
kf = []; m  = 6; n =18;
kf.m = m;  kf.n = n;
kf.Pk = 1*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db;[1;1;1]])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb; zeros(3,1); imuerr.sqg;imuerr.sqa; zeros(3,1);]).^2;
kf.Rk = diag(vperrset(1,2.5)).^2;
kf.Gammak = eye(kf.n); kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
kf.adaptive = 0; 
kf.pconstrain = 0;
kf.I = eye(kf.n);
%% others parameter setting
kf.adaptive = 1; 
kf.b = 0.5;kf.beta = 1;
kf.Rmin = diag(Rmin);
kf.Rmax = 100*kf.Rk;

kf.pconstrain = 1;
kf.Pkmax = (diag(kf.Pk)+1)*1e10; 
kf.Pkmin = Pmin;
ki =1;
gnss(find(gnss(:,7)>7),:) = [];
imugpssyn(imu(:,7), gnss(:,end));
%% Memory allocation
kf.Kk = zeros(kf.n,kf.m);
kf.Xk = zeros(kf.n,1);
kf.Phikk_1 = zeros(kf.n,kf.n);

t(1) = imu(1,end);
avp = zeros(N,10); avp(1,:) = [att',vel',pos',t(1)];
xkpk = zeros(N ,2*kf.n);   xkpk(1,:) = [kf.Xk',diag(kf.Pk)'];
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
        [kgps, ~] = imugpssyn(i, i, 'F');
        if kgps>0  && kgps <= L_GNSS
            kf.Xkk_1 = kf.Xk;
            kf.Pkk_1 = kf.Pk;
            Zk = [vel-gnss(kgps,1:3)';
                  pos-gnss(kgps,4:6)'];	
%             kf.Hk =zeros(kf.m,kf.n);kf.Hk(1:6,4:9) = eye(6);
			Mpv = [0,                 1/eth.RMh, 0;
				   1/(eth.RNh*eth.cl),0,         0;
				   0,                 0,         1];
			kf.Hk = [zeros(3,3) eye(3)     zeros(3,3) zeros(3,6) -Cnb*askew(wbib-Cnb'*eth.wnie);                                               % 
					 zeros(3,3) zeros(3,3) eye(3)     zeros(3,6) -Mpv*Cnb];
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
			if kf.pconstrain ==1
				for k = 1:kf.n
					if kf.Pk(k,k)<kf.Pkmin(k)
						kf.Pk(k,k) = kf.Pkmin(k);
					elseif kf.Pk(k,k)>kf.Pkmax
						ratio = sqrt(kf.Pkmax(k)/kf.Pk(k,k));
						kf.Pk(:,k) = kf.Pk(:,k)*ratio;
						kf.Pk(k,:) = kf.Pk(k,:)*ratio;
					end
				end
			end 			
            kf.Pk = (kf.Pk+kf.Pk')/2;
            [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua);             
            zk(ki,:) = [Zk',diag(kf.Rk)',ki];
            ki = ki+1;
        end
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
kf_res = varpack(avp,xkpk);
end