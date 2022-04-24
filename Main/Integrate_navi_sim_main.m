clc;clear all;
disp('抗差算法仿真验证');
%%  
load('trjRbs.mat');  
global glv 
glv = globalParameter(trjRbs.avp(1,7),trjRbs.avp(1,9));
%%  Change Dataformat to Standard
imu = trjRbs.imu;
L = length(imu);
t = imu(:,end);
gps = [trjRbs.gnss(:,4:6),trjRbs.gnss(:,1:3),repmat(2,length(trjRbs.gnss),1),trjRbs.gnss(:,7)];
avpmg = interp1(trjRbs.imu(:,end), trjRbs.avp,gps(:,8),'linear');  %线性插值
%添加误差
gps(320:332,1:6) =  gps(320:332,1:6) +5*(avpmg(320:332,4:9)-gps(320:332,1:6));
gps(450:470,1:6) =  gps(450:470,1:6) +10*(avpmg(450:470,4:9)-gps(450:470,1:6));
errgps = fplot([avpmg(:,1:3) gps(:,[1:6,end])],[avpmg,gps(:,8)],0);
%%  正常组合导航
avp0 = trjRbs.avp0;
davp0 = avpseterr([-2160;2520;-600], [0.3;0.3;0.3],[2;2;2]);%用来设置P,速度位置的参数对kf是最佳的，上面参数就不行
imuerr = imuerrset(3, 2000, trjRbs.arw, trjRbs.vrw);

%%Algorithm development
tic
resArkf = RobustKf(imu,gps,davp0,imuerr,avp0);
toc

tic
resBayes = BayesianFilter(imu,gps,davp0,imuerr,avp0);
toc

resMrkf = SH_KF156(imu,gps,davp0,imuerr,avp0);
resTkf = kf(imu,gps,davp0,imuerr,avp0);

%/-------------------------可以补充Filter_Algorithm中的其他算法用于对比------------------------------/


%%计算误差
errArkf  = fplot(resArkf .avp,[trjRbs.avp,trjRbs.imu(:,end)],0);
errBayes = fplot(resBayes.avp,[trjRbs.avp,trjRbs.imu(:,end)],0);
errTkf =  fplot(resTkf .avp,[trjRbs.avp,trjRbs.imu(:,end)],0);
errMrkf =  fplot(resMrkf .avp,[trjRbs.avp,trjRbs.imu(:,end)],0);

%%绘图比较
figure;
set(gcf,'position',[480 300 433 288]);
set(gca,'looseInset',[0 0 0 0]) %去空白区域
plot(t, trjRbs.avp(:,3)/glv.deg,'LineWidth',2);  axis([0 550 -200 200]);
ylabel(labeldef('y')); 
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
legend('Yaw');

figure
set(gcf,'position',[480 300 433 288]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
set(gca,'looseInset',[0 0 0 0]) %去空白区域
plot(t, [trjRbs.avp(:,4:6),sqrt(trjRbs.avp(:,4).^2+trjRbs.avp(:,5).^2+trjRbs.avp(:,6).^2)],'LineWidth',2); axis([0 550 -10 10]);
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel(labeldef('V')); h4=legend('V_E','V_N', 'V_U', '|V|');
set(h4,'Orientation','horizon')

figure
set(gcf,'position',[480 300 433 288]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
set(gca,'looseInset',[0 0 0 0]) %去空白区域
dxyz = pos2dxyz(trjRbs.avp(:,7:9));
plot(0, 0, 'rp'); 
hold on
plot(dxyz(:,1), dxyz(:,2),'LineWidth',2); 
xlabel(labeldef('est'));
ylabel(labeldef('nth'));
legend(sprintf('LON0:%.2f, LAT0:%.2f (DMS)', r2dms(trjRbs.avp(1,8)),r2dms(trjRbs.avp(1,7))));

%% 
figure('Color',[1 1 1]);
set(gcf,'position',[480 300 450 320]);
subplot(211),plot(errgps.errvel,'Linewidth',2);axis([0 550 -5 5]); % 设置坐标轴在指定的区间%set(gca,'XTick',[0:50:550]);
ylabel('\fontsize{10}\fontname{Times New Roman}Velocity Error(m/s)')
h4=legend('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_E',...
 '\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_N',...
  '\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_U');
set(h4,'Orientation','horizon')
subplot(212),plot(errgps.errpos,'Linewidth',2);axis([0 550 -50 50]);%set(gca,'XTick',0:50:550);
ylabel('\fontsize{10}\fontname{Times New Roman}Position Error (m)')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') ;%fontsize用来设置字体大小，fontname用来设置字体
h3=legend('\fontsize{10}\fontname{Times New Roman}\delta\itL\rm',...
 '\fontsize{10}\fontname{Times New Roman}\delta\it\lambda\rm',...
  '\fontsize{10}\fontname{Times New Roman}\delta\itH\rm');
set(h3,'Orientation','horizon')
%% 
figure('Color',[1 1 1]);
set(gcf,'position',[480 300 450 450]);
subplot(311),h1 = plot(t,errTkf.erratt(:,1),t,errMrkf .erratt(:,1),t,errBayes.erratt(:,1),t,errArkf.erratt(:,1),'LineWidth',2);axis([0 550 -0.8 0.2]);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_E/(\circ)')
subplot(312),h2 =plot(t,errTkf.erratt(:,2),t,errMrkf .erratt(:,2),t,errBayes.erratt(:,2),t,errArkf.erratt(:,2),'LineWidth',2);axis([0 550 -0.5 1]);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_N/(\circ)')
subplot(313),h3 =plot(t,errTkf.erratt(:,3),t,errMrkf .erratt(:,3),t,errBayes.erratt(:,3),t,errArkf.erratt(:,3),'LineWidth',2);axis([0 550 -20 2]);
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\it\phi\rm_U/(\circ)')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
h4 =legend('\fontsize{10}\fontname{Times New Roman}SKF',...
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
subplot(311),h1 = plot(t,errTkf.errpos(:,1),t,errMrkf .errpos(:,1),t,errBayes.errpos(:,1),t,errArkf.errpos(:,1),'LineWidth',2);axis([0 550 -5 8]);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itL\rm/ m')
subplot(312),h2 =plot(t,errTkf.errpos(:,2),t,errMrkf .errpos(:,2),t,errBayes.errpos(:,2),t,errArkf.errpos(:,2),'LineWidth',2);axis([0 550 -3 4]);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\it\lambda\rm/ m')
subplot(313),h3 =plot(t,errTkf.errpos(:,3),t,errMrkf .errpos(:,3),t,errBayes.errpos(:,3),t,errArkf.errpos(:,3),'LineWidth',2);axis([0 550 -2 2]);
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
set(gcf,'position',[480 300 450 450]);
subplot(311),h1 =plot(t,errTkf.errvel(:,1),t,errMrkf .errvel(:,1),t,errBayes.errvel(:,1),t,errArkf.errvel(:,1),'LineWidth',2);axis([0 550 -0.2 0.4]);
h1_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_E/ ( m/s )')
subplot(312),h2 =plot(t,errTkf.errvel(:,2),t,errMrkf .errvel(:,2),t,errBayes.errvel(:,2),t,errArkf.errvel(:,2),'LineWidth',2);axis([0 550 -0.4 0.4]);
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_N/ ( m/s )')
subplot(313),h3 =plot(t,errTkf.errvel(:,3),t,errMrkf .errvel(:,3),t,errBayes.errvel(:,3),t,errArkf.errvel(:,3),'LineWidth',2);axis([0 550 -0.1 0.2]);
h3_axis = gca; 
h4=legend('\fontsize{10}\fontname{Times New Roman}SKF',...
    '\fontsize{10}\fontname{Times New Roman}AKF',...
    '\fontsize{10}\fontname{Times New Roman}BayesFilter',...
    '\fontsize{10}\fontname{Times New Roman}ARKF');
ylabel('\fontsize{10}\fontname{Times New Roman}\delta\itV\rm_U/ ( m/s )')
xlabel('\fontsize{10}\fontname{Times New Roman}Time(s)'); %fontsize用来设置字体大小，fontname用来设置字体
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)
set(h4,'Orientation','horizon')

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
