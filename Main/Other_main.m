% this main is tested other function
clc;
clear all;
load('RealData.mat');
%%  Change Dataformat to Standard
imu = RealData.imu;
gps = RealData.gps;
gps = [gps(:,1:6) repmat(2,length(gps),1),gps(:,end)];
trj = RealData.trj;
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
t = imu(:,end);
%% 
avp0 = trj.avp0;
imuerr = imuerrset(3, 2000, 0.1, 100);
davp0 = avpseterr([60;60;60], [0.1;0.1;0.1], [10;10;10]); %初始误差
%% 
resNHC = NHC(imu,davp0,imuerr,avp0);
res  =  UKF156(imu,gps,davp0,imuerr,avp0);
res_all  = smooth_RTS_all(imu,gps,davp0,imuerr,avp0);
RTS_res = smooth_RTS(imu,gps,davp0,imuerr,avp0);
TKF_res = smooth_TKF(imu,gps,davp0,imuerr,avp0);
%% 

% err=fplot(res.avp,trj.avp,1); %均以弧度单位传入