function  err = fplot(avp,refavp,flag)
%% -----------Introduction------------
%通用画图函数
%input: 
%-------avp : 以弧度传入 N*10
%-------refavp : 以弧度传入 N*10
global glv
t = refavp(:,end);
set(0,'defaultfigurecolor','w') %figure 背景白色
set(gcf,'unit','centimeters','position',[13 9 12 9]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
%参考轨迹
% figure('name','avp参考运动信息')
% subplot(321), plot(t, refavp(:,1:2)/glv.deg,'Linewidth',2); xygo('pr');%航向和俯仰姿态角
% subplot(322), plot(t, refavp(:,3)/glv.deg,'Linewidth',2); xygo('y');%横滚姿态角
% subplot(323), plot(t, [refavp(:,4:6),sqrt(refavp(:,4).^2+refavp(:,5).^2+refavp(:,6).^2)],'Linewidth',2); xygo('V');%速度及合速度
% dxyz = dposxyz(refavp(:,7:9));
% subplot(325), plot(t, dxyz(:,1:3),'Linewidth',2); xygo('DP');%位置增量
% subplot(3,2,[4,6]), plot(0, 0, 'rp');   % 起点
% hold on, plot(dxyz(end,2), dxyz(end,1),'gp');   %终点
% hold on, plot(dxyz(:,2), dxyz(:,1),'Linewidth',2); xygo('est', 'nth');

%误差图
avp = [avp(:,1:3)/glv.deg avp(:,4:6) avp(:,7:8)/glv.deg avp(:,9) avp(:,10)];%化为标准单位
refavp = [refavp(:,1:3)./glv.deg refavp(:,4:6) refavp(:,7:8)./glv.deg refavp(:,9) refavp(:,10)];
erratt = -aa2phi(avp(:,1:3).*glv.deg,refavp(:,1:3).*glv.deg)./glv.deg;
errvel = avp(:,4:6) - refavp(:,4:6);
errpos = avp(:,7:9) - refavp(:,7:9);
errpos(:,1:2) = errpos(:,1:2).*60.*glv.nm;

err = varpack(erratt,errvel,errpos);
if flag
figure('name', 'KF姿态误差');
set(gcf,'unit','centimeters','position',[13 9 12 9]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
plot(t,erratt(:,1:3),'Linewidth',2);grid on;
xlabel('\fontsize{12}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel('\fontsize{12}\fontname{Times New Roman}Attitude Error (deg)')
legend('\fontsize{12}\fontname{宋体}俯仰角','\fontsize{12}\fontname{宋体}横滚角','\fontsize{12}\fontname{宋体}航偏角');

figure('name', 'KF速度误差');
set(gcf,'unit','centimeters','position',[13 9 12 9]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
plot(t,errvel(:,1:3),'Linewidth',2);grid on;
xlabel('\fontsize{12}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel('\fontsize{12}\fontname{Times New Roman}Velocity Error (m/s)')
legend('\fontsize{12}\fontname{宋体}东向速度','\fontsize{12}\fontname{宋体}北向速度','\fontsize{12}\fontname{宋体}天向速度');

figure('name', 'KF位置误差');
set(gcf,'unit','centimeters','position',[13 9 12 9]);
set(gcf, 'Color', [1,1,1]);%图标外围设为白色
plot(t,errpos(:,1:3),'Linewidth',2);grid on;
xlabel('\fontsize{12}\fontname{Times New Roman}Time(s)') %fontsize用来设置字体大小，fontname用来设置字体
ylabel('\fontsize{12}\fontname{Times New Roman}Position Error (m)')
legend('\fontsize{12}\fontname{宋体}纬度','\fontsize{12}\fontname{宋体}经度','\fontsize{12}\fontname{宋体}高度');
end
end