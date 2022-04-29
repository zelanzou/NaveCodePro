%% Algorithm Introduction
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author :
% Version :
% Date : 
% File : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all;                                                                 % close all figures
clc;                                                                       % clear cmd text
clear ;                                                                    % clear all RAM
disp('杆臂估计算法验证');
%% Load data
load('trj.mat');   
glvs
%% Algorithm time
m = length(trj.imu);                                                       % 数据量
n = (1:m)';                                                                % 数据标号
dt = trj.ts;                                                               % 采样时间
t = trj.imu(:,7);  % 数据时长
imuerr = imuerrset(1, 500, 0.01, 40);
% imuerr = imuerrset(0.03, 100, 0.001, 10);
imu = imuadderr(trj.imu, imuerr);
imu(:,1:6) = imu(:,1:6)/dt;
% 2 将杆臂误差添加在gps中
davp0 = avpseterr([30;-30;30], [0.1;0.1;0.1], [1;1;1]);
lever = [1.5; -0.5; 1.0];
gps = gpssimu(trj.avp, davp0(4:6), davp0(7:9), 1, lever);%在GPS仿真数据中加入杆臂误差
%gps = gpssimu(trj.avp, davp0(4:6), davp0(7:9)); 不添加杆臂

%验证imu输出
figure('name','陀螺仪输出')
plot(t,imu(:,1),t,imu(:,2),t,imu(:,3));
figure('name','加速度计输出')
plot(t,imu(:,4),t,imu(:,5),t,imu(:,6));

%% Algorithm develop

% res186 = Sins_Gps_186(imu,gps,davp0,imuerr,trj);
res186 = Virtual_Lever_Arm(imu,gps,davp0,imuerr,trj);
%% Error analysis
%% Plot figures
%% Save Workspace
%% End