clc;
clear ;
load('test9250_04.txt')
dt=1/200;
pos0 = [31.303996*pi/180;120.6195853*pi/180;12.34];  %实验室当地经纬度
global glv 
glv = globalParameter(pos0(1),pos0(3));
imu =[ test9250_04(:,2:4), test9250_04(:,5:7)*glv.g0,test9250_04(:,1)];
%% 标加计

zup = imu(7*200+1:12*200,1:6);
xup = imu(61*200+1:66*200,1:6);
zdown = imu(128*200+1:133*200,1:6);
xdown = imu(145*200+1:149*200,1:6);
yup = imu(216*200+1:221*200,1:6);
ydown =imu(249*200+1:254*200,1:6);
vs = [xup; xdown; yup; ydown;zup; zdown];
accVector = [mean(xup(:,4:6));mean(xdown(:,4:6));mean(yup(:,4:6));mean(ydown(:,4:6));mean(zup(:,4:6));mean(zdown(:,4:6))]';

resAcc0301 = mAcce6PosCalibration(accVector(:,1),accVector(:,2),accVector(:,3),accVector(:,4),accVector(:,5),accVector(:,6),glv.g0);

resAcc0301.AsbL(1,2)=0;  resAcc0301.AsbL(1,3)=0; resAcc0301.AsbL(2,3)=0; %构造下三角模型e

%% 标陀螺零偏

gybias = mean(vs(:,1:3));
imu(:,1:3) = imu(:,1:3)  - gybias;

%2.2 标比例因子
zup_plusR = imu(30*200+1:40*200,1:6); %20秒
zup_minusR = imu(42*200+1:52*200,1:6); %20秒
gzup = ( trapz(zup_plusR)- trapz(zup_minusR))/200;

xdown_plusR= imu(65*200+1:77*200,1:6); %10秒
xdown_minusR  = imu(78*200+1:90*200,1:6); %10秒
gxdown = (trapz(xdown_plusR) - trapz(xdown_minusR ) )/200;

zdown_plusR = imu(118*200+1:130*200,1:6); %10秒
zdown_minusR = imu(105*200+1:117*200,1:6); %10秒
gzdown =  (trapz(zdown_plusR) - trapz(zdown_minusR ) )/200;

xup_plusR = imu(148*200+1:158*200,1:6); %10秒
xup_minusR = imu(158*200+1:168*200,1:6); %10秒
gxup= (trapz(xup_plusR) - trapz(xup_minusR ))/200;

yup_plusR = imu(202*200+1:212*200,1:6); %12秒
yup_minusR = imu(192*200+1:202*200,1:6); %12秒
gyup= ( trapz(yup_plusR) - trapz(yup_minusR  ))/200;

ydown_plusR =  imu(250*200+1:270*200,1:6); %10秒
ydown_minusR = imu(233*200+1:253*200,1:6);%10秒
gydown = ( trapz(ydown_plusR) - trapz(ydown_minusR ))/200;

Rot = [zup_plusR;zup_minusR;xup_plusR ;xup_minusR;zdown_plusR ;zdown_minusR;...
xdown_plusR;xdown_minusR ;yup_plusR ;yup_minusR;ydown_plusR ;ydown_minusR];

gyVector = [gxup(1:3);gxdown(1:3);gyup(1:3);gydown(1:3);gzup(1:3);gzdown(1:3)]';

resGyro=mGyro6PosCalibration(gyVector(:,1),gyVector(:,2),gyVector(:,3),gyVector(:,4),gyVector(:,5),gyVector(:,6),720);

%2.3 标非正交误差，算法对比

data = imu(285*200:end,:);


%已标定参数数据补偿

fb = (resAcc0301.SFm*resAcc0301.AsbL)^-1 *(data(:,4:6) - resAcc0301.bisv')';
gbl = resGyro.SFm^-1*data(:,1:3)';
% gbl = resGyro.SFm^-1*(data(:,1:3)-gybias)';

pitch0 = asin(mean(fb(2,1:800))/glv.g0);
roll0 = atan2(-mean(fb(1,1:800)), mean(fb(3,1:800)));
avp0=[pitch0;roll0;0;zeros(3,1);pos0];
rer=Arw_Vrw2std(0.5, 225, dt,pos0,1) ;      % 消费级 MEMS imu参数
kfpara.ab=[1;1;1]*rer.acc_std;      % 加计静态标准差
kfpara.wb=[1;1;1]*rer.gyro_std*pi/180;     % 陀螺静态标准差
davp = avpseterr([720;720;30], [0.01;0.01;0.01], [0.05;0.05;0.05]); %初始误差
imuerr.web =kfpara.wb;  % 0.02*glv.dpsh;
imuerr.wdb = kfpara.ab;   %100*glv.ugpsHz;
kf.Pk=diag([0.2*pi/180*ones(1,3) 0.05*ones(1,3) 0.05*ones(1,3)]).^2;
kf.Qk=diag(1.5*kfpara.wb(1)*ones(1,3)).^2;
kf.Rk=diag(1.5*kfpara.wb(1)*glv.g0*ones(1,3));
kf.Xk=zeros(9,1);

%算法1
CalibrateResult1 =kfCalibrated_Xu([gbl'*pi/180 fb' (1/200:1/200:length(fb)/200)'],[],imuerr,avp0);
Gsb = [1,CalibrateResult1.xkpk(end,10:11); ...
            CalibrateResult1.xkpk(end,12),1,CalibrateResult1.xkpk(end,13);...
            CalibrateResult1.xkpk(end,14:15),1]    ;
 res1=Gsb*(resGyro.SFm^-1);
%算法2

CalibrateResult2 = mfunRef43(data(:,1:3)'*pi/180, fb, (1/200:1/200:length(fb)/200)',[],kf,pos0,2);

res2 = eye(3,3)+diag([CalibrateResult2.Xki(end,4),CalibrateResult2.Xki(end,5),CalibrateResult2.Xki(end,6)])-...
                        askew([CalibrateResult2.Xki(end,7),CalibrateResult2.Xki(end,8),CalibrateResult2.Xki(end,9)]); 
%% step3 验证标定结果
%测试数据
load('Data\test9250.txt')

imu0 = [test9250(1000:13000,2:4),test9250(1000:13000,5:7)*glv.g0,(1/200:1/200:length(test9250(1000:13000,:))/200)'];
 %校正
fbib = (resAcc0301.SFm*resAcc0301.AsbL)^-1 *(imu0(:,4:6) - resAcc0301.bisv')';
wbib1 = (res1)^-1*(imu0(:,1:3) - gybias)';
wbib2 = (res2)^-1*(imu0(:,1:3) - gybias)';
 
pitch0 = asin(mean(fbib(2,1:3500)/glv.g0));
roll0 = atan2(-mean(fbib(1,1:3500)), mean(fbib(3,1:3500)));

avp0=[pitch0;roll0;0;zeros(3,1);pos0];
res1inpure = ins_inpure([ wbib1' *pi/180, fbib',imu0(:,end)],avp0);
res2inpure = ins_inpure([ wbib2' *pi/180, fbib',imu0(:,end)],avp0);
resInpurebefore =  ins_inpure([imu0(:,1:3)*pi/180, imu0(:,4:6),imu0(:,end)],avp0);

imuerr = imuerrset(3,1000,0.03,100);
gps = [repmat([0;0;0;pos0;2]',length(imu0),1),imu0(:,end)];
ref =  smooth_TKF([ imu0(:,1:3) *pi/180, imu0(:,4:end)],gps,davp,imuerr,avp0); %此处应该用已校正的数据

%% 
t = imu0(:,end);
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(t,ref.xs(:,1)*180/pi,'r','Linewidth',2);
hold on 
plot(t,res2inpure.avp(:,1)*180/pi,'b:','Linewidth',2);grid on
hold on 
plot(t,res1inpure.avp(:,1)*180/pi,'m--','Linewidth',2);grid on
legend('\fontsize{10}\fontname{宋体}参考值','\fontsize{10}\fontname{宋体}传统方法','\fontsize{10}\fontname{宋体}两步修正法');
ylabel('\fontsize{10}\fontname{宋体}俯仰角(°)')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(t,ref.xs(:,2)*180/pi,'r','Linewidth',2);
hold on 
plot(t,res2inpure.avp(:,2)*180/pi,'b:','Linewidth',2);grid on
hold on 
plot(t,res1inpure.avp(:,2)*180/pi,'m--','Linewidth',2);grid on
legend('\fontsize{10}\fontname{宋体}参考值','\fontsize{10}\fontname{宋体}传统方法','\fontsize{10}\fontname{宋体}两步修正法');
ylabel('\fontsize{10}\fontname{宋体}横滚角(°)')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体


%% 
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(t,ref.xs(:,1)*180/pi,'r','Linewidth',2);
hold on 
plot(t,resInpurebefore.avp(:,1)*180/pi,'b:','Linewidth',2);grid on
hold on 
plot(t,res1inpure.avp(:,1)*180/pi,'m--','Linewidth',2);grid on
legend('\fontsize{10}\fontname{宋体}参考值','\fontsize{10}\fontname{宋体}标定前','\fontsize{10}\fontname{宋体}标定后');
ylabel('\fontsize{10}\fontname{宋体}俯仰角(°)')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(t,ref.xs(:,2)*180/pi,'r','Linewidth',2);
hold on 
plot(t,resInpurebefore.avp(:,2)*180/pi,'b:','Linewidth',2);grid on
hold on 
plot(t,res1inpure.avp(:,2)*180/pi,'m--','Linewidth',2);grid on
legend('\fontsize{10}\fontname{宋体}参考值','\fontsize{10}\fontname{宋体}标定前','\fontsize{10}\fontname{宋体}标定后');
ylabel('\fontsize{10}\fontname{宋体}横滚角(°)')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体
%% 
marker_idx = 1:5000:length(CalibrateResult1.avp(:,end));
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,10),'bx-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,11),'g+-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,12),'rs-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,13),'cd-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,14),'mo-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,15),'yp-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
grid on;
legend('\fontsize{10}\fontname{Times New Roman}k12','\fontsize{10}\fontname{Times New Roman}k13','\fontsize{10}\fontname{Times New Roman}k21','\fontsize{10}\fontname{Times New Roman}k23','\fontsize{10}\fontname{Times New Roman}k31','\fontsize{10}\fontname{Times New Roman}k32')
ylabel('\fontsize{10}\fontname{宋体}非正交误差估计值')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,25),'bx-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,26),'g+-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,27),'rs-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,28),'cd-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,29),'mo-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(CalibrateResult1.avp(:,end),CalibrateResult1.xkpk(:,30),'yp-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
grid on;
ylabel('\fontsize{10}\fontname{宋体}单位：度(°)')
legend('k12','k13','k21','k23','k31','k32')
legend({'$$\sqrt{P11}$$','$$\sqrt{P22}$$','$$\sqrt{P33}$$','$$\sqrt{P44}$$','$$\sqrt{P55}$$','$$\sqrt{P66}$$'},'interpreter','latex')
xlabel('\fontsize{10}\fontname{宋体}时间/s') %fontsize用来设置字体大小，fontname用来设置字体
%% 
marker_idx = 1:5000:length(data);
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(data(:,4),'bx-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(data(:,5),'go-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(data(:,6),'r--','MarkerIndices',marker_idx,'Linewidth',2);hold on;
grid on;
ylabel('\fontsize{10}\fontname{Times New Roman}m/s^2')
legend('\fontsize{10}\fontname{Times New Roman}x','\fontsize{10}\fontname{Times New Roman}y','\fontsize{10}\fontname{Times New Roman}z')
xlabel('\fontsize{10}\fontname{宋体}采样时刻')

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 10 7]);
set(gca,'looseInset',[0 0 0 0]) %去白边
plot(data(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(data(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);hold on;
plot(data(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);hold on;
grid on;
ylabel('\fontsize{10}\fontname{Times New Roman} °/s')
legend('\fontsize{10}\fontname{Times New Roman}x','\fontsize{10}\fontname{Times New Roman}y','\fontsize{10}\fontname{Times New Roman}z')
xlabel('\fontsize{10}\fontname{宋体}采样时刻')
