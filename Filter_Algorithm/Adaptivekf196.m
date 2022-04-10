function kf_res = Adaptivekf196(imu,gnss,davp,imuerr,avp0,Rmin,Pmin)
%% -----------Introduction------------
%19维自适应卡尔曼滤波,考虑同步时间
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
kf = []; m  = 6; n =19;
kf.m = m;  kf.n = n;
kf.Pk = 1*diag([davp(1:3); davp(4:6); davp(7:9); imuerr.eb; imuerr.db;[1;1;1];0.1])^2;
kf.Qt = diag([imuerr.web; imuerr.wdb; zeros(3,1); imuerr.sqg;imuerr.sqa; zeros(3,1);0]).^2;
kf.Rk = diag(vperrset(1,2.5)).^2;
kf.Gammak = eye(kf.n); kf.Gammak(1:3,1:3) = -Cnb; kf.Gammak(4:6,4:6) = Cnb;
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
timebar(1,N,'19维量测自适应+方差受限');
for i= 2:N 
        wbib = imu(i,1:3)' ;
        fb = imu(i,4:6)' ;
        dt = imu(i,end) - imu(i-1,end);
        t(i) = imu(i,end);
%         an = Cnb*fb + eth.gcc;
        %% 惯导更新
%         [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
            eth = EarthParameter(pos,vel);
            fn = Cnb*fb;
            an = fn + eth.gcc; %同步时间的加速度
            vel = vel + an*dt;

            Mpv = [0 1/eth.RMh 0;1/eth.clRNh 0 0;0 0 1];
            Mpvvn = Mpv*vel;
            pos = pos + Mpvvn*dt;

            eth = EarthParameter(pos,vel);
            wbnb = wbib - Cnb'*eth.wnin;%b系下
            qua = qupdt(qua,wbnb*dt);%wbnb*dt是相当于等效旋转增量rv
            Cnb = q2mat(qua);
            att = m2att(Cnb);
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
% 			Mpv = [0,                 1/eth.RMh, 0;
% 				   1/(eth.RNh*eth.cl),0,         0;
% 				   0,                 0,         1];
%             Mpvvn  = Mpv*vel;
			kf.Hk = [zeros(3,3) eye(3)     zeros(3,3) zeros(3,6) -Cnb*askew(wbib-Cnb'*eth.wnie) -an;                                               % 
					 zeros(3,3) zeros(3,3) eye(3)     zeros(3,6) -Mpv*Cnb                       -Mpvvn];
            kf.rk = Zk - kf.Hk*kf.Xkk_1;%残差
            kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差
			kf.Py0 =  kf.Hk*kf.PXZkk_1;
			if kf.adaptive==1  % for adaptive KF, make sure Rk is diag 24/04/2015 ？?
				for k=1:kf.m
					ry = kf.rk(k)^2 -  kf.Py0(k,k); %序贯量测+方差受限
					if ry<kf.Rmin(k,k),  kf.Rk(k,k) = kf.Rmin(k,k); end
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
            [att,pos,vel,qua,Cnb,kf,xk]= feed_back_correct(kf,[att;vel;pos],qua); 
%             kf.xk(10:15) = 0;
            zk(ki,:) = [Zk',diag(kf.Rk)',ki];
            ki = ki+1;
        end
        avp(i,:) = [att',vel',pos',t(i)];
        xkpk(i,:) = [xk',diag(kf.Pk)'];
        timebar;
end
kf_res = varpack(avp,xkpk,zk);
end