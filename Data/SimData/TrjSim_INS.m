%% Algorithm Introduction
% The row vector is recommended
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The trajectory simulation script.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author :
% Version : V0.1.0
% Date : 
% File : 
% output:
%		trj.ts 采样时间                      
%		trj.avp 导航系下无误差的导航参数 
%			单位：rad  rad  rad  m/s^2   rad rad m s 
%		 		[pitch roll yaw ve vn vu lat lon h t]
%		trj.imu  惯性单元输出
%				 角增量   速度增量 s 
%				[wx,wy,wz,ax,ay,az,t]
%		trj.avp0 初始导航参数
%		trj.wat  运动轨迹参数
%            	wat(:,1) - period lasting time 持续时间
%            	wat(:,2) - period initial velocity magnitude 初始速度
%            	wat(:,3:5) - angular rate in trajectory-frame 角速率
%            	wat(:,6:8) - acceleration in trajectory-frame 加速度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;                                                                       % clear command
close all;                                                                 % close figures
clear ;                                                                    % clear data space
disp('This is a trajectory simulation script');
%% Initial data 
glvs
ts = 0.005;       %  200HZ sampling interval 总时长 722s
avp0 = avpset([0;0;0], [0,0,0], glv.pos0); % init avp. if need the initial velocity, it is need to set in avp0
%初始速度最好是为0 ，否则含有除速度做转弯轨迹会有问题。
%% Function index
xxx = [];
%% Algorithm develop
% trajectory segment setting

% seg = trjsegment(xxx, 'init',         0);
% seg = trjsegment(seg, 'uniform',      50);
% seg = trjsegment(seg, 'accelerate',   20, xxx, 1.5);%以1m/s^2加速 
% seg = trjsegment(seg, 'uniform',      100);
% seg = trjsegment(seg, 'turnleft',   45, 2, xxx, 4);%协调左转弯时，
% seg = trjsegment(seg, 'uniform',      135);
% seg = trjsegment(seg, 'deaccelerate',   5, xxx, 1);%以1m/s^2加速 
% % seg = trjsegment(seg, 'coturnright',  45, 2, xxx, 4);%协调右转弯，
% seg = trjsegment(seg, 'turnright',   45, 2, xxx, 4);%协调左转弯时，
% seg = trjsegment(seg, 'uniform',      100);
% seg = trjsegment(seg, 'deaccelerate', 5, xxx, 1);%以1m/s^2减速 
% seg = trjsegment(seg, 'uniform',      100);
% seg = trjsegment(seg, 'turnright',   45, 2, xxx, 4);%协调左转弯时，
% seg = trjsegment(seg, 'uniform',      100);
% seg = trjsegment(seg, 'accelerate', 5, xxx, 1);%以1m/s^2减速 
% seg = trjsegment(seg, 'uniform',      100);
% seg = trjsegment(seg, 'turnleft',   45, 2, xxx, 4);%协调左转弯时，
% seg = trjsegment(seg, 'uniform',      100);
% % seg = trjsegment(seg, 'deaccelerate',   5, xxx, 1);%以1m/s^2加速 
% trj = trjsimu(avp0, seg.wat, ts, 1);                 % truth data

% seg = trjsegment(xxx, 'init',         0);
% seg = trjsegment(seg, 'uniform',      25*60);
% trj = trjsimu(avp0, seg.wat, ts, 1);
%% 
% 8字运动
% seg = trjsegment(xxx, 'init',         0);
% seg = trjsegment(seg, 'uniform',     60);
% seg = trjsegment(seg, 'accelerate',   5, xxx, 1); 
% seg = trjsegment(seg, '8turn',        90, 2, xxx, 4);
% seg = trjsegment(seg, 'uniform',      10);
% seg = trjsegment(seg, 'deaccelerate', 5,  xxx, 1);
% seg = trjsegment(seg, 'uniform',      1500);
% trj = trjsimu(avp0, seg.wat, ts, 1);                 % truth data

%% 
%椭圆运动
% seg = trjsegment(xxx, 'init',         0);%必须先对seg初始化一次
% seg = trjsegment(seg, 'uniform',     10);
% seg = trjsegment(seg, 'coturnleft',   45, 2,xxx, 4);
% seg = trjsegment(seg, 'coturnleft',   45, 2,xxx, 4);
% seg = trjsegment(seg, 'coturnleft',   45, 2,xxx, 4);
% seg = trjsegment(seg, 'coturnleft',   45, 2,xxx, 4);
% trj = trjsimu(avp0, seg.wat, ts, 1);                 % truth data

%% 
%%动静态交替估计
% xxx = [];
% seg = trjsegment(xxx, 'init',         0);
% seg = trjsegment(seg, 'uniform',      10);
% seg = trjsegment(seg, 'accelerate',   4, xxx, 5);%以1m/s^2加速 
% seg = trjsegment(seg, 'uniform',      20);
% seg = trjsegment(seg, 'deaccelerate',   5, xxx, 4);%以1m/s^2加速 
% seg = trjsegment(seg, 'uniform',      10);
% seg = trjsegment(seg, 'accelerate',   5, xxx, 7);%以1m/s^2加速 
% seg = trjsegment(seg, 'deaccelerate',   10, xxx, 3);%以1m/s^2加速 
% seg = trjsegment(seg, 'uniform',      10);
% seg = trjsegment(seg, 'deaccelerate',   5, xxx, 1);%以1m/s^2加速 
% seg = trjsegment(seg, 'uniform',      20);
% seg = trjsegment(seg, 'accelerate',   10, xxx, 1);%以1m/s^2加速 
% trj = trjsimu(avp0, seg.wat, ts, 1);                 % truth data

%% 
%转弯模式数据
% xxx = [];
seg = trjsegment(xxx, 'init',         0);
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'accelerate',   4, xxx, 5);%以1m/s^2加速 
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'turnleft',   3, 2, xxx, 4);%左转弯时，小转弯
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'deaccelerate', 5, xxx, 2);%以1m/s^2减速 
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'turnright',   10, 9, xxx, 4);%大右转弯
seg = trjsegment(seg, 'uniform',      20);
seg = trjsegment(seg, 'accelerate',   3, xxx, 2);%以1m/s^2加速
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'turnleft',   60, 3, xxx, 4);%左转弯时，小转弯
seg = trjsegment(seg, 'uniform',      30);
seg = trjsegment(seg, 'turnright',   2, 5, xxx, 4);%小右转弯
seg = trjsegment(seg, 'uniform',      30);
trj = trjsimu(avp0, seg.wat, ts, 1);                 % truth data

%% Plot figures
insplot(trj.avp);
imuplot(trj.imu);
%% Save Workspace
% pos(rad,rad,m) vel(m/s) att(rad)
save('trj.mat', 'trj');


%% End