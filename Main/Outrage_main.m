clc
clear all;
disp('GNSS Outrage算法验证');
%% Load data
load('imuSTGnss0608_F1.mat')
%%  Change Dataformat to Standard
index = 1:length(imuSTGn0608_1.ref(:,1));
gindex = 1:length(imuSTGn0608_1.gnss(:,1));
load('f1savp.mat') %补偿安装角之后保存的子惯导数据，可做参考值
trj.avp = savp; %补偿了主子惯导安装角，可以作为参考值
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
imu = [imuSTGn0608_1.imu1(index,2:4)*glv.deg imuSTGn0608_1.imu1(index,5:7)*glv.g0 imuSTGn0608_1.imu1(index,12)];
gps = [imuSTGn0608_1.gnss(gindex,7:9) imuSTGn0608_1.gnss(gindex,4:5)*pi/180 imuSTGn0608_1.gnss(gindex,6) imuSTGn0608_1.gnss(gindex,11) imuSTGn0608_1.gnss(gindex,10)];
gps = [interp1(gps(:,end), gps(:,1:end-1),[0:1:2829]','linear'),[0:1:2829]'];  %线性插值 原始数据漏了122和1972两帧

index = [610,660,1501,1680,2100,2250];
gps([610:660,1501:1680,2100:2250],:) = [];%模拟中断情况

t = imu(:,end);
avp0 = trj.avp(1,:)';

imuerr = imuerrset(3, 1000, 0.03, 100);
davp0 = avpseterr([60;-60;60], [1;1;1], [2.5;2.5;2.5]); %初始误差
Rmin
Rmax
%-----------------------可以加入Outrage_Solution中的其他算法---------------------/
result = VmP_NHC_RTS(imu,gps,davp,imuerr,avp0);

errres = fplot(res.avp,trj.avp,1); %均以弧度单位传入