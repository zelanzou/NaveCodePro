function splot(bf,af,navs)
%% 画图 
figure(1);
plot(bf,'r','Linewidth',1);
hold on;
plot(af,'b','Linewidth',1);
legend('最佳适应值','平均适应值')
xlabel('种群代数');ylabel('适应度值');
figure(2);
subplot(3,1,1);
plot(navs(:,1),'r','Linewidth',1);
hold on;
plot(navs(:,2),'g','Linewidth',1);
hold on;
plot(navs(:,3),'b','Linewidth',1);
legend('Bx','By','Bz');
xlabel('种群代数');ylabel('零偏值/(m*s-2)');
% figure(3);
subplot(3,1,2)
plot(navs(:,4),'r','Linewidth',1);
hold on;
plot(navs(:,5),'g','Linewidth',1);
hold on;
plot(navs(:,6),'b','Linewidth',1);
legend('Sxx','Syy','Szz');
xlabel('种群代数');ylabel('标度因数值');
% figure(4);
subplot(3,1,3)
plot(navs(:,7),'r','Linewidth',1);
hold on;
plot(navs(:,8),'g','Linewidth',1);
hold on;
plot(navs(:,9),'b','Linewidth',1);
legend('Mxx','Myy','Mzz');
xlabel('种群代数');ylabel('安装误差值');
 end