function res = ins_inpure(imu,avp0)
%% -----------Introduction------------
%纯惯性导航
%input: 
%-------imu : 传感器数据N*7 单位：rad/s   m/s^2
%-------avp0  ； 初始姿态信息  9*1
%output
%-------res.avp: N*10 导航信息，标准单位弧度、m/s 、m 
%%  data length and time
N = length(imu(:,end));
%% Memory allocation                                                    
avp = zeros(N,10);
%% Initial data
att = avp0(1:3);vel = avp0(4:6);pos = avp0(7:9);
qua = a2qua(att);
Cnb = a2mat(att);
t(1) = imu(1,end);
avp(1,:) = [att',vel',pos',t(1)];
timebar(1,N,'ins inpure navigation');
for i = 2:N
	wbib = imu(i,1:3)' ;
	fb = imu(i,4:6)' ;
	dt = imu(i,end) - imu(i-1,end);
	t(i) = imu(i,end);
    [att,vel,pos,qua,Cnb] = avp_update(wbib,fb,Cnb,qua,pos,vel,dt);
	avp(i,:) = [att',vel',pos',t(i)];
	timebar;
end
%% Save Workspace
res =  varpack(avp);
end