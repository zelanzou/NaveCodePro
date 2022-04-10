function res = UKF_constrait(imu,davp,imuerr,trj)
%% Function Introduction
% UKF + 内部约束：包含高程，速度，水平姿态角
%% Function input

%% Function time
%% Algorithm time
dt = imu(2,7)-imu(1,7);
L_INS = length(imu);
t = imu(:,end);
nts = 1* dt;
%% Memory allocation                                                    
avp = zeros(L_INS,10);
xkpk = zeros(L_INS,30);
%% Initial data 
glvs
att = trj.avp0(1:3);vel = trj.avp0(4:6);pos = trj.avp0(7:9);
avp(1,:) = [att',vel',pos',t(1)];
Vu = vel(3);h = pos(3);pitch = att(1);roll = att(2);vb0  = [0;0;0];
qua = a2qua(att);
Cnb = a2mat(att);

kf = [];
kf.Pk = 1*diag([1*davp(1:3); 1*davp(4:6); 1*davp(7:9); 1*imuerr.eb; 1*imuerr.db])^2;
kf.Qk = diag([ imuerr.web; imuerr.wdb;zeros(9,1)])^2*dt;
kf.Rk = 10*diag([0.05;0.1;0.05;0.05;0.2;0.1])^2;%   
kf.Hk = zeros(6,15); 
[kf.m, kf.n] = size(kf.Hk);
kf.Kk = zeros(kf.n, kf.m);
kf.Xk = zeros(kf.n, 1);
kf.Gammak = eye(kf.n);
kf.I = eye(kf.n);
xkpk(1,:) = [kf.Xk',diag(kf.Pk)']; 
%UKF参数初始化
[gamma,Wm,Wc] = UKFParameter(kf.n);
ki = timebar(1, L_INS, 'UKF+内部约束模型验证.');

for i = 2:L_INS
	wbib = imu(i,1:3)';
	fb = imu(i,4:6)';
	[att,vel,pos,eth,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
	Vu = [Vu vel(3)];h = [h,pos(3)];pitch = [pitch,att(1)];roll = [roll,att(2)];
	Cbn = Cnb';
	vb = Cbn*vel;vb0 = [vb0,vb];
	M2 = -Cbn*askew(vel);M1 = Cbn;	
	Z = [
		  Vu(i)-Vu(i-1);
		  h(i)-h(i-1);
		  vb(1);
		  vb(3);
		  pitch(i)-pitch(i-1);
		  roll(i)-roll(i-1);
		];
	H = [
		 zeros(1,5) 1       zeros(1,9);
		 zeros(1,8) 1       zeros(1,6);
		 M2(1,:)    M1(1,:) zeros(1,9);
		 M2(3,:)    M1(3,:) zeros(1,9);
		 1          0       zeros(1,13);
		 0          1       zeros(1,13);
		];
    kf = ukf_filter(kf,[att;vel;pos],imu,gamma,Wm,Wc,nts,1,Z,H);
	xk = kf.Xk;
    %进行反馈
    vel = vel - kf.Xk(4:6);
    pos = pos- kf.Xk(7:9);
    qua = qdelphi(qua, kf.Xk(1:3));%qdelafa 是大失准角
    Cnb = q2mat(qua);att = m2att(Cnb);
    kf.Xk(1:9)= 0;   
	xkpk(i,:) = [xk',diag(kf.Pk)']; 
    avp(i,:) = [att',vel',pos',t(i)];
    timebar;
end 		
%% Error analysis
% avp(:,4:9) = avp(:,4:9) - xkpk(:,4:9);
avp = [avp(:,1:3)/glv.deg avp(:,4:6) avp(:,7:8)/glv.deg avp(:,9) avp(:,10)];%化为标准单位
trj.avp = [trj.avp(:,1:3)./glv.deg trj.avp(:,4:6) trj.avp(:,7:8)./glv.deg trj.avp(:,9) trj.avp(:,10)];	
erratt = -aa2phi(avp(:,1:3)*glv.deg,trj.avp(:,1:3)*glv.deg)/glv.deg; %计算失准角，应为小角
errvel = avp(:,4:6) - trj.avp(:,4:6);
errpos = avp(:,7:9) - trj.avp(:,7:9);
errpos(:,1:2) = errpos(:,1:2).*60.*glv.nm;
%% Plot figures
set(0,'defaultfigurecolor','w') %figure 背景白色
figure('name', 'UKF姿态误差');
plot(t,erratt(:,1),'r',t,erratt(:,2),'g',t,erratt(:,3),'b');grid on;
figure('name', 'UKF速度误差');
plot(t,errvel(:,1),'r',t,errvel(:,2),'g',t,errvel(:,3),'b');grid on;
figure('name', 'UKF位置误差');
plot(t,errpos(:,1),'r',t,errpos(:,2),'g',t,errpos(:,3),'b');grid on;
%% Save Workspace
res =  varpack(avp,xkpk); 
end
