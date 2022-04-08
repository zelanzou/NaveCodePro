function res=mfunRef43(wib_b,fib_b,tm,imupara,kf,pos0,flag)
% 文献中的转动顺序及时间. 该仿真数据比较理想,没有实际手动转动时的干扰、振动.
%--------------------------------------------------------------------------
% 输入参数：
%               imupara:结构体{Asba, Asbg, Sa, Sg, abias, gbias, wb,ab}
%               Asba, Asbg: 正交系到加计、陀螺的变换矩阵.     Sa,Sg: 加计、陀螺标度因数
%               abias,gbias: 加计、陀螺零偏, 3*1, 单位：m/s2, rad/s
%               ab,wb:  加计、陀螺静态数据标准差. 可根据随机游走和函数 Arw_Vrw2std()获取.
%  Cb0p: 固定imu的框架,放置在大理石平台上, 框架初始时,加计下三角正交系 与平台系之间变换矩阵(包含加计与框架之间的安装角)
%  g0: 当地 重力数值.      L0:当地纬度       dt:采样间隔
%  例如: imupara.Asba=eye(3,3);  imupara.Asbg=eye(3,3)+askew([0.5;0.5;0.6]*pi/180);
%        imupara.Sa=eye(3,3)*0.9992;  imupara.Sg=eye(3,3)*0.98; 
% imupara.abias=[1;1;1]*0.03; imupara.gbias=[1;1;1]*0.05;   imupara.ab=[1;1;1]*0.01; imupara.wb=[1;1;1]*0.05;
%  输出参数：
%    结构体res=varpack(Xki,GyroBias,AccBias,posi,atti,Veni,stPki,gmInd);
%  
% Author : xutongxu
% Version : 1 
% Date :2021.6.27
% File: Zhou Q ,  Yu G ,  Li H , et al. A Novel MEMS Gyroscope In-Self Calibration Approach[J]. Sensors, 2020, 20(18).
%% 初始化参数
% 至少前9s数据为静止
% eth=earth(pos0,zeros(3,1));
eth = EarthParameter(pos0,zeros(3,1));
dt=tm(2)-tm(1);   nt=round(1/dt);
fb0=mean(fib_b(:,1:nt*4),2);
wb0=mean(wib_b(:,1:nt*4),2);       % 静止陀螺零偏

pitch0=asin(fb0(2)/eth.g);
roll0=atan2(-fb0(1), fb0(3));

att0=[pitch0; roll0; 0];  Cnb=a2mat(att0);  qnb=m2qua(Cnb);
N=length(tm);

% --------KF 参数初始化---------
if(flag==2)
    % 外部给出 kf参数
    kf.Xk=zeros(9,1);
else
    % 自动给出
    kf.Pk=diag([3*imupara.ab(1)*ones(1,3)/eth.g, ones(1,3)*0.06, ones(1,3)*0.08]).^2;
    kf.Qk=diag(1.5*imupara.wb(1)*ones(1,3)).^2;
    kf.Rk=diag(1.5*imupara.wb(1)*eth.g*ones(1,3));
    kf.Xk=zeros(9,1);
end
Hk=zeros(3,9);   Hk(2,1)=eth.g;    Hk(1,2)=-eth.g;    G=zeros(9,3);
%% 分配内存
Xki=zeros(N,9);   
Pki=zeros(N,9);   Zki=zeros(3,N);   atti=zeros(3,N);
eIn=eye(9,9);     atti(:,1)=att0;
%% 文献内算法
g0=eth.g;
gn=([0,0,-g0])'; 
fn= -gn;
% 由于前面数据静止,所以从第2帧开始估计. 参照文献, 数据为摇摆数据集, 几乎无线加速度.

for i=2:N
     dt=tm(i)-tm(i-1);
     wbibk=0.5*(wib_b(:,i)+wib_b(:,i-1))-wb0;
     fib=fib_b(:,i);
     % 1. 姿态更新
     qnb=qupdt(qnb,wbibk*dt);
     Cnbk=q2mat(qnb);        
   
     % 2. 计算量测
     fnj=(Cnbk)*fib;
     Zk=fnj-fn;
     Zki(:,i)=Zk;
     % 3. KF滤波
     Ft=getFt(Cnb, wbibk);  G(1:3,1:3)=Cnb;  Gk=G*dt;
     Phik=eIn +Ft*dt;
     
     Xkk_1=Phik*kf.Xk;     
     Pkk_1=Phik*kf.Pk*(Phik') + Gk*kf.Qk*(Gk');
%      Pkk_1 =  0.5* (Pkk_1 + Pkk_1');       
     Kk=Pkk_1*(Hk')*((Hk*Pkk_1*(Hk') + kf.Rk)^(-1));
     % 估计
     kf.Xk=Xkk_1 + Kk*(Zk - Hk*Xkk_1);
     kf.Pk=(eIn - Kk*Hk)*Pkk_1;
     
     % 4. 反馈
     % 姿态、航向修正
%      ra=1;
%      qx=rv2q(ra*kf.Xk(1:3));    Cn_n1=q2mat(qx);
%      Cnbk=(Cn_n1)*Cnbk;         
%      kf.Xk(1:3)=kf.Xk(1:3)*(1-ra);
    
     Xki(i,:)=kf.Xk';
     Pki(i,:)=diag(kf.Pk)';
     
     Cnb=Cnbk;  qnb=m2qua(Cnb);   atti(:,i)=m2att(Cnb);
end
Xki(1,:)=Xki(2,:);
Pki(1,:)=Pki(2,:);
res=varpack(Xki,Pki,Zki,atti);

end

function Ft=getFt(Cnb,wib_b)
Ft=zeros(9,9);
Ft(1:3,4:6)=Cnb*diag(wib_b);
Ft(1:3,7:9)=-Cnb*askew(wib_b);

end
