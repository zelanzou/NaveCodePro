function [imu,vsb,fb] = SimuAcc(phi,flag)
%% Function Introduction
% 24位置加速度计传感器输出数据仿真。
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variables:
%   The input is:
%		phi:水平倾角，用弧度表示
%       f : 绘图标志位（1：绘图，2：不绘图）
%   The output as:
%       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : zouzelan
% Version : V2.0 
% Date : 2020.5.29
% File : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
glvs;
%% Function input

%% Function time
dt=0.0025;                                            % Sampling time  400HZ
tf=5;                                                % Simulation time
t=dt:dt:tf;                                          % Time stamp
n = length(t);
pos = [31.8897 118.82033 12.34];%苏州地区的经、纬、高度
pos = pos.*(pi/180);%化为弧度
eth = earth(pos);
gn = -eth.gn;
%% Memory allocation
imu = [];fb_cross = [];
%% Initial data 
%设真值
Sa = diag([0.99, 0.99102 ,0.99024]);
bias =[0.14736; 0.14736; 0.14736];
Csb = [1      0       0
       pi/180  1       0
       pi/180 pi/180   1];%非正交误差
Cph = [ 1  0        0
        0  cos(phi) -sin(phi)
        0  sin(phi) cos(phi)];
Chn = eye(3,3);
switch(flag)
    case 6
        P = [0,-pi,0,0,-pi/2,pi/2];
        R = [0,0,pi/2,-pi/2,0,0];
        H = [0,0,0,0,0,0];
    case 24
        P = [0 0 0 0 -pi -pi -pi -pi 0 pi/2 pi -pi/2 0 pi/2 pi -pi/2 -pi/2 -pi/2 -pi/2 -pi/2 pi/2 pi/2 pi/2 pi/2];
        R = [0 0 0 0 0 0 0 0 pi/2 pi/2 pi/2 pi/2 -pi/2 -pi/2 -pi/2 -pi/2 0 pi/2 pi -pi/2 0 pi/2 pi -pi/2];
        H = [0 pi/2 pi -pi/2 0 pi/2 pi -pi/2 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
    case 8
        P = [0 0 0 0 0 0 0 0];
        R = [0 0 0 0 0 0 0 0];
        H = [0 pi/4 pi/2 pi*3/4 pi pi*5/4 pi*3/2 pi*7/4];
end
% ki = timebar(1, 24, 'imu Simulation Algorithm.');                                   % 第一个参数表示间隔数，第二个参数表示总进度

%% Function index

%% Function realize
for i = 1:length(P)    
    noise = 0.005*randn(n,3);%24个位置取不同的噪声
    beta = (rand(1)-0.5)*2;%消除手动转动的误差,垂直转动无误差，水平转动有误差
    cbeta = cos(beta);sbeta = sin(beta);
    Czbeta = [cbeta sbeta 0;sbeta cbeta 0;0 0 1];
    Cxbeta = [1 0 0;0 cbeta sbeta;0 sbeta cbeta];
    Cybeta = [cbeta 0 sbeta;0 1 0;sbeta 0 cbeta];
    cp=cos(P(i));sp=sin(P(i));cr=cos(R(i));sr=sin(R(i));ch=cos(H(i));sh=sin(H(i));
    Cb2 = [cr  0  sr
           0   1  0
           -sr 0  cr];%绕y轴转
    C21 = [1 0  0 
           0 cp -sp
           0 sp cp];%绕x轴转
    C1p = [ch -sh  0
           sh ch  0
           0  0   1];%绕z轴转 
    if flag == 6
        if i == 1||i == 2|| i == 5 || i == 6
            Cbp = C1p*C21*Cb2;
        else
            Cbp = Cb2*C21*C1p;
        end
    elseif flag == 24
            if i==1 || i==5 || i==9 || i==13
                Cbp = C1p*C21*Cb2;
            end
            if i>=2&&i<=8
                Cbp = Czbeta*C1p*C21*Cb2;
            end
            if i>=10&&i<=16
                 Cbp = Cxbeta*C21*Cb2*C1p;
            end
            if i==17 || i==21
               Cbp = Cb2*C21*C1p;
            end
            if i>17&&i<=24
               Cbp = Cybeta*Cb2*C21*C1p;
            end
    elseif flag == 8
        Cbp =  Czbeta*Cb2*C21*C1p;
    end
    Csn = Csb*Cbp*Cph*Chn;
    fs= Csn*(gn);%传感器系下的重力矢量
    fb = Sa * (fs - bias);%正交系下的重力矢量
    fs = repmat(fs',n,1) + noise;
    fb = repmat(fb',n,1)+noise; 
    imu = [imu; fs];
    fb_cross = [fb_cross;fb];
%     ki = timebar;                                                          % 进度条
end

%计算均值
 for j = 1:length(P) 
    vsb(j,:) = mean(imu(2000*(j-1)+1:2000*j,:));
 end
%  vsb = vsb';
%% Error Analysis

%% Plot figures
plot(imu);
legend('x轴','y轴','z轴');
%% Return Results and Save Workspace
save imu.mat imu
save vsb.mat vsb


end

%% End

