clc;clear all;
%% Load data
load('trj.mat');  
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
%%  Change Dataformat to Standard
imuerr = imuerrset(3, 3000, 0.02, 70); %Ìí¼ÓÆ÷¼þÎó²î
imu = imuadderr(trj.imu, imuerr);
imu(:,1:6) = imu(:,1:6)/trj.ts;
lever = [0; 0; 0];
gps = gpssimu(trj.avp, [0.2;0.2;0.2], posseterr([2;2;2]));% ²»Ìí¼Ó¸Ë±Û
gps = [gps(:,1:6) repmat(2,length(gps),1) gps(:,7)];
%% Initial data
davp0 = avpseterr([60;-60;60], [0.2;0.2;0.2], [2;2;2]); %³õÊ¼Îó²î
imuerr = imuerrset(3, 3000, 0.02, 70);
avp0 = avpadderr(trj.avp0,davp0);
Rmin = vperrset(0.1, 0.1).^2;
Pmin = [avperrset([0.2,3],0.01,0.01);gabias(0.1, [30,30]); [0.01;0.01;0.01]; 0.001].^2;
%% Algorithm develop
resAdaptivekf196 =  Adaptivekf196(imu,gps,davp0,imuerr,avp0,Rmin,Pmin);
err=fplot(resAdaptivekf196.avp,trj.avp,1);