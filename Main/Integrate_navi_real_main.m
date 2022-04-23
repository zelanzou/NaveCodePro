%% Algorithm Introduction
clc;
clear all;
disp('自适应抗差算法跑车数据验证');
%% Load data
load('imuSTGnss0608_F1.mat')
%%  Change Dataformat to Standard
index = 1:length(imuSTGn0608_1.ref(:,1));
gindex = 1:length(imuSTGn0608_1.gnss(:,1));
load('E:\跑车数据\BeijingGnssIMUata\BeijingCAR\f1savp.mat') %补偿安装角之后保存的子惯导数据，可做参考值
trj.avp = savp; %补偿了主子惯导安装角，可以作为参考值
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
imu = [imuSTGn0608_1.imu1(index,2:4)*glv.deg imuSTGn0608_1.imu1(index,5:7)*glv.g0 imuSTGn0608_1.imu1(index,12)];
gps = [imuSTGn0608_1.gnss(gindex,7:9) imuSTGn0608_1.gnss(gindex,4:5)*pi/180 imuSTGn0608_1.gnss(gindex,6) imuSTGn0608_1.gnss(gindex,11) imuSTGn0608_1.gnss(gindex,10)];
t = imu(:,end);
avp0 = trj.avp(1,:)';

imuerr = imuerrset(3, 1000, 0.03, 100);
davp0 = avpseterr([60;-60;60], [0.1;0.1;0.1], [2.5;2.5;2.5]); %初始误差
%%  添加粗差
%方式，用已有的gps误差，等倍数放大
avpmg = interp1(trj.avp(:,end), trj.avp,gps(:,8),'linear');  %线性插值

errgps = fplot([avpmg(:,1:3) gps(:,[1:6,end])],[avpmg,gps(:,8)],1);

gps(1000:1050,1:6) = gps(1000:1050,1:6) +5*(avpmg(1000:1050,4:9)- gps(1000:1050,1:6));

gps(2200:2220,1:6) = gps(2200:2220,1:6)  +8*(avpmg(2200:2220,4:9)- gps(2200:2220,1:6));

gps(1500:1510,6) =gps(1500:1510,6) +10*(avpmg(1500:1510,6)- gps(1500:1510,6));

%%Algorithm development
tic
resArkf = RobustKf(imu,gps,davp0,imuerr,avp0);
toc

tic
resBayes = BayesianFilter(imu,gps,davp0,imuerr,avp0);
toc
%resMrkf = ARkf(imu,gps,davp0,imuerr,avp0);
resMrkf = SH_KF156(imu,gps,davp0,imuerr,avp0);
resTkf =  kf(imu,gps,davp0,imuerr,avp0);

%/-------------------------可以补充Filter_Algorithm中的其他算法用于对比------------------------------/


%%计算误差
errArkf  = fplot(resArkf .avp,trj.avp,1);
errBayes = fplot(resBayes.avp,trj.avp,0);
errTkf =  fplot(resTkf .avp,trj.avp,0);
errMrkf =  fplot(resMrkf .avp,trj.avp,0);

%%绘图比较
figure('Color',[1 1 1]);
set(gcf,'position',[480 300 450 450]);
subplot(311),h1 = plot(t,errTkf.erratt(:,1),t,errMrkf .erratt(:,1),t,errBayes.erratt(:,1),t,errArkf.erratt(:,1),'LineWidth',2);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_E/(\circ)')
subplot(312),h2 =plot(t,errTkf.erratt(:,2),t,errMrkf .erratt(:,2),t,errBayes.erratt(:,2),t,errArkf.erratt(:,2),'LineWidth',2);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_N/(\circ)')
subplot(313),h3 =plot(t,errTkf.erratt(:,3),t,errMrkf .erratt(:,3),t,errBayes.erratt(:,3),t,errArkf.erratt(:,3),'LineWidth',2);
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_U/(\circ)')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF');
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)
% set(h4,'Orientation','horizon')

figure('Color',[1 1 1]);
% set(gcf,'unit','centimeters','position',[0 0 8.76 6.64]);
set(gcf,'position',[480 300 450 450]);
subplot(311),h1 = plot(t,errTkf.errpos(:,1),t,errMrkf .errpos(:,1),t,errBayes.errpos(:,1),t,errArkf.errpos(:,1),'LineWidth',2);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itL\rm/ m')
subplot(312),h2 =plot(t,errTkf.errpos(:,2),t,errMrkf .errpos(:,2),t,errBayes.errpos(:,2),t,errArkf.errpos(:,2),'LineWidth',2);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\it\lambda\rm/ m')
subplot(313),h3 =plot(t,errTkf.errpos(:,3),t,errMrkf .errpos(:,3),t,errBayes.errpos(:,3),t,errArkf.errpos(:,3),'LineWidth',2);
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itH\rm/ m')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
h4=legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF');
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)
set(h4,'Orientation','horizon')

figure('Color',[1 1 1]);
% set(gcf,'unit','centimeters','position',[0 0 8.76 6.64]);
set(gcf,'position',[480 300 450 450]);
subplot(311),h1 =plot(t,errTkf.errvel(:,1),t,errMrkf .errvel(:,1),t,errBayes.errvel(:,1),t,errArkf.errvel(:,1),'LineWidth',2);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_E/ ( m/s )')
subplot(312),h2 =plot(t,errTkf.errvel(:,2),t,errMrkf .errvel(:,2),t,errBayes.errvel(:,2),t,errArkf.errvel(:,2),'LineWidth',2);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_N/ ( m/s )')
subplot(313),h3 =plot(t,errTkf.errvel(:,3),t,errMrkf .errvel(:,3),t,errBayes.errvel(:,3),t,errArkf.errvel(:,3),'LineWidth',2);
h3_axis = gca; 
legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF')
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_U/ ( m/s )')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)

%% 
figure;
set(gcf,'position',[480 300 433 288]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
set(gca,'looseInset',[0 0 0 0]) %去空白区域
plot(t, trj.avp(:,3)/glv.deg,'LineWidth',2);  
ylabel(labeldef('y')); 
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
legend('Yaw');

figure
set(gcf,'position',[480 300 433 288]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
set(gca,'looseInset',[0 0 0 0]) %去空白区域
plot(t, [trj.avp(:,4:6),sqrt(trj.avp(:,4).^2+trj.avp(:,5).^2+trj.avp(:,6).^2)],'LineWidth',2); 
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel(labeldef('V')); 
h4=legend('V_E','V_N', 'V_U', '|V|');
set(h4,'Orientation','horizon')

%% 
figure;
% set(gcf,'unit','centimeters','position',[0 0 10.83 8.21]);
set(gcf,'position',[480 300 433 288]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
set(gca,'looseInset',[0 0 0 0]) %去空白区域
plot(gps(:,7),'LineWidth',2);
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel('\fontsize{10}\fontname{Times New Roman}PDOP')

%% 
figure('Color',[1 1 1]);
set(gcf,'position',[480 300 450 320]);
subplot(211),h1 =plot(errgps.errvel,'Linewidth',2); 
 h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}Velocity Error(m/s)')
h4=legend('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_E',...
 '\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_N',...
  '\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_U');
set(h4,'Orientation','horizon')
subplot(212),h2 =plot(errgps.errpos,'Linewidth',2);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}Position Error (m)')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') ;
legend('\fontsize{10}\fontname{Times New Roman}\delta\itL\rm',...
 '\fontsize{10}\fontname{Times New Roman}\delta\it\lambda\rm',...
  '\fontsize{10}\fontname{Times New Roman}\delta\itH\rm');
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)

%计算RMSE
[RMSETkf,maxerrTkf] = rmse(errTkf.errpos);
[RMSEBayes,maxerrBayes] = rmse(errBayes.errpos);
[RMSEArkf,maxArkf] = rmse(errArkf .errpos);
[RMSEMrkf,maxMrkf] = rmse(errMrkf .errpos);
%% bar 绘图
A=[RMSETkf;RMSEMrkf;RMSEBayes;RMSEArkf]';
B = [maxerrTkf;maxMrkf;maxerrBayes;maxArkf]';

figure('name','RMSE')
set(0,'defaultfigurecolor','w') %figure 背景白色
set(gcf,'position',[480 300 433 288]);
bar(A)
legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF')
ylabel('\fontsize{10}\fontname{Times New Roman}RMSE(m)')
set(gca,'XTickLabel',{'\fontsize{10}\fontname{Times New Roman}Latitude' '\fontsize{10}\fontname{Times New Roman}Longtiude' '\fontsize{10}\fontname{Times New Roman}Altitude'});

figure('name','最大误差')
set(0,'defaultfigurecolor','w') %figure 背景白色
set(gcf,'position',[480 300 433 288]);
bar(B)
legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF')
ylabel('\fontsize{12}\fontname{Times New Roman}Maximum Error Absolute  (m)')
set(gca,'XTickLabel',{'\fontsize{10}\fontname{Times New Roman}Latitude' '\fontsize{10}\fontname{Times New Roman}Longtiude' '\fontsize{10}\fontname{Times New Roman}Altitude'});
