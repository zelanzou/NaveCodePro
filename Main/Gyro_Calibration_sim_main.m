clear;clc; 
glvs
Re=6378137.0;                   f=0.003352813177897;
%% 
pos0 = [31.303996*pi/180;120.6195853*pi/180;12.34];
dt=1/2000;
imupara.Asba=eye(3,3) ;  % 可改为MEMS 加计能标定到的精度
imupara.Asbg=eye(3,3) ;
imupara.Sa=eye(3,3);                             % 可改为 MEMS 加计可标定到的精度,0.9998  
imupara.Sg=eye(3,3) ; 
imupara.abias=[0.0;0.0;0.0];           % 可修改为加计标定后剩余误差 0.0015
imupara.gbias=[0.0;0.0;0.0];           % 零偏 假设陀螺零偏已经校正，可修改为还剩于不可标定误差
imupara.ab=[0;0;0];      % 加计静态标准差
imupara.wb=[0;0;0];     % 陀螺静态标准差
eth=earth(pos0,zeros(3,1));      g0=eth.g;
r=[0.05;0.05;0.05];
Cb0b = a2mat([0.3;0.4;0.5]*pi/180);
%Ng=-askew([0.025; -0.02; 0.01]);  
eS=0.15*[0.05; -0.04; -0.03];
%Asbg=eye(3,3) + Ng;
ang1 = 0.45;ang2 = 0.65;
accAngle = [0,0,-ang1,0,-ang1,ang1]*pi/80;
gyroAngle = [-ang1,ang2,-ang1,ang2,-ang1,ang2]*pi/80;
[Asba,Asbg] = genImuAsb(gyroAngle,accAngle);  % 以摇摆方式产生矩阵
Sg=eye(3,3) +diag(eS)*1; 
gbias=[0.5;0.5;0.5]*1*pi/180; 
%% 
rer=Arw_Vrw2std(0.5, 225, dt,pos0,1) ;      % 消费级 MEMS imu参数
kfpara.ab=[1;1;1]*rer.acc_std;      % 加计静态标准差
kfpara.wb=[1;1;1]*rer.gyro_std*pi/180;     % 陀螺静态标准差
davp = avpseterr([720;720;30], [0.01;0.01;0.01], [0.05;0.05;0.05]); %初始误差
imuerr.web =kfpara.wb;  % 0.02*glv.dpsh;
imuerr.wdb = kfpara.ab;   %100*glv.ugpsHz;

kf.Pk=diag([0.2*pi/180*ones(1,3) 0.05*ones(1,3) 0.05*ones(1,3)]).^2;
kf.Qk=diag(1.5*kfpara.wb(1)*ones(1,3)).^2;
kf.Rk=diag(1.5*kfpara.wb(1)*eth.g*ones(1,3));

%% 
imuData = cell(100,1);
res1= cell(100,1);
res2 = cell(100,1);
k1mat = zeros(100,9);
k2mat = zeros(100,9);
%% 
for i=1:100
    att = [-8;-8;3]+i*[0.08;0.08;0]; %范围在-8~8°范围线性变化
    Cb0p=a2mat(att*pi/180)'; %放置imu与框架之间的夹角。
    rsDa=genGyro_Data(imupara,Cb0p,g0,pos0(1),dt);
    tm=(1:length(rsDa.twsbib(1,:)))*dt;
    for j=10:length(tm)
        dwbib = (rsDa.twsbib(:,j)-rsDa.twsbib(:,j-1))/dt;
        rsDa.tfsbib(:,j)=Cb0b * (rsDa.tfsbib(:,j) +cross(dwbib,r)+cross(rsDa.twsbib(:,j),cross(rsDa.twsbib(:,j),r)));  %假设存在向心加速度
		rsDa.twsbib(:,j) = Cb0b*rsDa.twsbib(:,j);
    end
    fsb = rsDa.tfsbib(:,11:end);  tm(:,1:10) = [];
    wsb = rsDa.twsbib(:,11:end);
    imuData{i} = adderr(fsb,wsb,tm,Asbg,Sg,gbias);

    
    imu = imuData{i} ;
    resGyro = scaleCali(imu);
    
    data = imu(402*100+1:end,:);
    fb =data(:,4:6)';
    gbl = resGyro.SFm^-1*data(:,1:3)';
    pitch0 = asin(mean(fb(2,1:400))/g0);
    roll0 = atan2(-mean(fb(1,1:400)), mean(fb(3,1:400)));
    avp0=[pitch0;roll0;0;zeros(3,1);pos0];
    
    CalibrateResult1 = kfCalibrated_Xu([gbl'*pi/180 fb' (1/100:1/100:length(fb)/100)'],davp,imuerr,avp0);
    CalibrateResult2 = mfunRef43(data(:,1:3)'*pi/180, fb, (1/100:1/100:length(fb)/100)',[],kf,pos0,2);
    
    Gsb = [1,CalibrateResult1.xkpk(end,10:11); ...
                CalibrateResult1.xkpk(end,12),1,CalibrateResult1.xkpk(end,13);...
                CalibrateResult1.xkpk(end,14:15),1]    ;
    res1{i}=Gsb*(resGyro.SFm^-1);
    res2{i}=eye(3,3)+diag([CalibrateResult2.Xki(end,4),CalibrateResult2.Xki(end,5),CalibrateResult2.Xki(end,6)])-...
                        askew([CalibrateResult2.Xki(end,7),CalibrateResult2.Xki(end,8),CalibrateResult2.Xki(end,9)]) ;
     k1mat(i,:) = [res1{i}(1,:), res1{i}(2,:),res1{i}(3,:)];
     k2mat(i,:)= [res2{i}(1,:), res2{i}(2,:),res2{i}(3,:)];    

end
ref = (Sg*Asbg)^(-1);

%% 画图结果分析
figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 15 11]);
a =  [k1mat(:,1),k2mat(:,1),repmat(ref(1,1),100,1)];
marker_idx = 1:15:length(a);
subplot(311); h1 =plot( a(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(a(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(a(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2) ;grid on
h1_axis = gca;
legend('\fontsize{10}\fontname{宋体}两步修正法','\fontsize{10}\fontname{宋体}传统方法','\fontsize{10}\fontname{宋体}参考真值');
ylabel('\fontsize{10}\fontname{Times New Roman}k11')
b= [k1mat(:,5),k2mat(:,5),repmat(ref(2,2),100,1)];
subplot(312); h2 =plot( b(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(b(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(b(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k22')
c= [k1mat(:,9),k2mat(:,9),repmat(ref(3,3),100,1)];
subplot(313); h3 =plot( c(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(c(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(c(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k33')
xlabel('\fontsize{10}\fontname{宋体}重复次数/次') %fontsize用来设置字体大小，fontname用来设置字体
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 15 11]);
a =[k1mat(:,2),k2mat(:,2),repmat(ref(1,2),100,1)];
subplot(311); h1 =plot( a(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(a(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(a(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2) ;grid on
h1_axis = gca; 
legend('\fontsize{10}\fontname{宋体}两步修正法','\fontsize{10}\fontname{宋体}传统方法','\fontsize{10}\fontname{宋体}参考真值');
ylabel('\fontsize{10}\fontname{Times New Roman}k12')
b = [k1mat(:,3),k2mat(:,3),repmat(ref(1,3),100,1)];
subplot(312); h2 =plot( b(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(b(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(b(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k13')
c = [k1mat(:,6),k2mat(:,6),repmat(ref(2,3),100,1)];
subplot(313); h3 =plot( c(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(c(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(c(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k23')
xlabel('\fontsize{10}\fontname{宋体}重复次数/次') %fontsize用来设置字体大小，fontname用来设置字体
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)

figure('Color',[1 1 1]);
set(gcf,'unit','centimeters','position',[13 9 15 11]);
a = [k1mat(:,4),k2mat(:,4),repmat(ref(2,1),100,1)];
subplot(311); h1 =plot( a(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(a(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(a(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2) ;grid on
h1_axis = gca; 
legend('\fontsize{10}\fontname{宋体}两步修正法','\fontsize{10}\fontname{宋体}传统方法','\fontsize{10}\fontname{宋体}参考真值');
ylabel('\fontsize{10}\fontname{Times New Roman}k21')
b=[k1mat(:,7),k2mat(:,7),repmat(ref(3,1),100,1)];
subplot(312); h2 =plot( b(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(b(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(b(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h2_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k31')
c = [k1mat(:,8),k2mat(:,8),repmat(ref(3,2),100,1)];
subplot(313); h3 =plot( c(:,1),'bx-','MarkerIndices',marker_idx,'Linewidth',2);
hold on; plot(c(:,2),'go-','MarkerIndices',marker_idx,'Linewidth',2);
hold on;plot(c(:,3),'r--','MarkerIndices',marker_idx,'Linewidth',2);grid on
h3_axis = gca; 
ylabel('\fontsize{10}\fontname{Times New Roman}k32')
xlabel('\fontsize{10}\fontname{宋体}重复次数/次') %fontsize用来设置字体大小，fontname用来设置字体
Expand_axis_fill_figure(h1_axis)
Expand_axis_fill_figure(h2_axis)
Expand_axis_fill_figure(h3_axis)

%%计算RMSE
k11 = sqrt(1/length(k1mat) *sum((k1mat(:,9)-ref(3,3)).^2));
sqrt(1/length(k2mat) *sum((k2mat(:,9)-ref(3,3)).^2));
%均值
mean1 = mean(k1mat- [ref(1,1),ref(1,2),ref(1,3),ref(2,1),ref(2,2),ref(2,3),ref(3,1),ref(3,2),ref(3,3)]);
mean2 = mean(k2mat- [ref(1,1),ref(1,2),ref(1,3),ref(2,1),ref(2,2),ref(2,3),ref(3,1),ref(3,2),ref(3,3)]);
%标准差
std1 =sqrt( sum((k1mat-mean(k1mat)).^2)/length(k1mat));
std2 =sqrt( sum((k2mat-mean(k2mat)).^2)/length(k2mat));
%一次相对误差
CalibrateResult1_err= [k1mat(1,1)-ref(1,1),k1mat(1,5)-ref(2,2),k1mat(1,9)-ref(3,3),...
						k1mat(1,2)-ref(1,2),k1mat(1,3)-ref(1,3),k1mat(1,6)-ref(2,3),...
						k1mat(1,4)-ref(2,1),k1mat(1,7)-ref(3,1),k1mat(1,8)-ref(3,2)];
CalibrateResult2_err= [k2mat(1,1)-ref(1,1),k2mat(1,5)-ref(2,2),k2mat(1,9)-ref(3,3),...
						k2mat(1,2)-ref(1,2),k2mat(1,3)-ref(1,3),k2mat(1,6)-ref(2,3),...
						k2mat(1,4)-ref(2,1),k2mat(1,7)-ref(3,1),k2mat(1,8)-ref(3,2)];
%% %相对误差
relative1 = (res1{1}- ref)*100./ref;
relative2 = (res2{1}- ref)*100./ref;              
%% 

function imu = adderr(fsb,wsb,tm,Asbg,Sg,gbias)
fsb = fsb+0.001+randn(3,length(tm))*0.015;
wsbib = Sg*Asbg*wsb +gbias + randn(3,length(tm))*0.05*pi/180;
imu = [wsbib(:,1:20:end)',fsb(:,1:20:end)',tm(1:20:end)'];
gybias = mean(imu(1:600,1:3));
imu(:,1:3) = imu(:,1:3)  - gybias    ;
imu(:,1:3) = imu(:,1:3) *180/pi;
end
%% 

function  resGyro = scaleCali(imu)
zup_plusR = imu(5*100+1:30*100,1:6); %20秒
zup_minusR = imu(30*100+1:55*100,1:6); %20秒
gzup = ( trapz(zup_plusR)- trapz(zup_minusR))/100;

xdown_plusR= imu(205*100+1:230*100,1:6); %10秒
xdown_minusR  = imu(230*100+1:255*100,1:6); %10秒
gxdown = (trapz(xdown_plusR) - trapz(xdown_minusR ) )/100;

zdown_plusR = imu(140*100+1:165*100,1:6); %10秒
zdown_minusR = imu(165*100+1:190*100,1:6); %10秒
gzdown =  (trapz(zdown_plusR) - trapz(zdown_minusR ) )/100;

xup_plusR = imu(75*100+1:100*100,1:6); %10秒
xup_minusR = imu(100*100+1:125*100,1:6); %10秒
gxup= (trapz(xup_plusR) - trapz(xup_minusR ))/100;

yup_plusR = imu(335*100+1:360*100,1:6); %12秒
yup_minusR = imu(360*100+1:385*100,1:6); %12秒
gyup= ( trapz(yup_plusR) - trapz(yup_minusR  ))/100;

ydown_plusR = imu(270*100+1:295*100,1:6); %10秒
ydown_minusR = imu(295*100+1:320*100,1:6); %10秒
gydown = ( trapz(ydown_plusR) - trapz(ydown_minusR ))/100;

Rot = [zup_plusR;zup_minusR;xup_plusR ;xup_minusR;zdown_plusR ;zdown_minusR;...
xdown_plusR;xdown_minusR ;yup_plusR ;yup_minusR;ydown_plusR ;ydown_minusR];

gyVector = [gxup(1:3);gxdown(1:3);gyup(1:3);gydown(1:3);gzup(1:3);gzdown(1:3)]';

resGyro=mGyro6PosCalibration(gyVector(:,1),gyVector(:,2),gyVector(:,3),gyVector(:,4),gyVector(:,5),gyVector(:,6),720);
end
