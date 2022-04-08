function [Asba,Asbg]= genImuAsb(gyroAngle,acceAngle)
% 产生非正交变换矩阵，从理想 b 系到 传感器系（陀螺仪、加速度计坐标系）
% 按右手转动顺序产生矩阵，例如 Zs: 等于Zb 先绕 x 轴转 -> 再绕 y轴 转
%   vector_s=Asba*vector_b; 
%   vector_s=Asbg*vector_b;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 输入变量：
%        gyroAngle:[xy,xz,yz,yx,zx,zy]   % 按右手转动顺序的角度，单位：弧度
%        acc3Angle:[xy,xz,yz,yx,zx,zy]
% 输出变量：
%        Asba: 加速度计正交系到非正交系  变换矩阵
%        Asbg: 陀螺仪正交系到非正交系  变换矩阵
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : xutongxu
% Version : 1
% Date : 2019.12.16
% File : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Csx_b=mgenRotM(gyroAngle(2),'z',2)*mgenRotM(gyroAngle(1),'y',2);
Rx=Csx_b(1,:);                               % b系矢量在 s 系 x轴上投影
Csy_b=mgenRotM(gyroAngle(4),'x',2)*mgenRotM(gyroAngle(3),'z',2);
Ry=Csy_b(2,:);                               % b系矢量在 s 系 x轴上投影
Csz_b=mgenRotM(gyroAngle(6),'y',2)*mgenRotM(gyroAngle(5),'x',2);
Rz=Csz_b(3,:);                               % b系矢量在 s 系 x轴上投影
Asbg=[Rx;Ry;Rz];                              % 正交系b系矢量，在s系投影 变换矩阵
Rxb=(Csx_b')*[1;0;0];
Ryb=(Csy_b')*[0;1;0];
Rzb=(Csz_b')*[0;0;1];
ax2y=acos(sum(Rxb.*Ryb)/(norm(Rxb)*norm(Ryb)))*180/pi;
ax2z=acos(sum(Rxb.*Rzb)/(norm(Rxb)*norm(Rzb)))*180/pi;
ay2z=acos(sum(Ryb.*Rzb)/(norm(Ryb)*norm(Rzb)))*180/pi;
fprintf('\n陀螺仪  ：x-y夹角：%.3f度, x-z夹角：%.3f度,y-z夹角：%.3f度\n',ax2y,ax2z,ay2z);

Csx_b=mgenRotM(acceAngle(2),'z',2)*mgenRotM(acceAngle(1),'y',2);
Rx=Csx_b(1,:);                               % b系矢量在 s 系 x轴上投影
Csy_b=mgenRotM(acceAngle(4),'x',2)*mgenRotM(acceAngle(3),'z',2);
Ry=Csy_b(2,:);                               % b系矢量在 s 系 x轴上投影
Csz_b=mgenRotM(acceAngle(6),'y',2)*mgenRotM(acceAngle(5),'x',2);
Rz=Csz_b(3,:);                               % b系矢量在 s 系 x轴上投影
Rxb=(Csx_b')*[1;0;0];
Ryb=(Csy_b')*[0;1;0];
Rzb=(Csz_b')*[0;0;1];
Asba=[Rx;Ry;Rz];                              % 正交系b系矢量，在s系投影 变换矩阵
ax2y=acos(sum(Rxb.*Ryb)/(norm(Rxb)*norm(Ryb)))*180/pi;
ax2z=acos(sum(Rxb.*Rzb)/(norm(Rxb)*norm(Rzb)))*180/pi;
ay2z=acos(sum(Ryb.*Rzb)/(norm(Ryb)*norm(Rzb)))*180/pi;
fprintf('加速度计：x-y夹角：%.3f度, x-z夹角：%.3f度,y-z夹角：%.3f度\n',ax2y,ax2z,ay2z);
% 
% Ka=Asba;
% Kg=Asbg;
% ax2y=acos(sum(Ka(:,1).*Ka(:,2))/(norm(Ka(:,1))*norm(Ka(:,2))))*180/pi;
% ax2z=acos(sum(Ka(:,1).*Ka(:,3))/(norm(Ka(:,1))*norm(Ka(:,3))))*180/pi;
% ay2z=acos(sum(Ka(:,2).*Ka(:,3))/(norm(Ka(:,2))*norm(Ka(:,3))))*180/pi;
% 
% fprintf('陀螺仪  ：x-y夹角：%.3f度, x-z夹角：%.3f度,y-z夹角：%.3f度\n',ax2y,ax2z,ay2z);
 end

function mrot=mgenRotM(radian,xyz,check)
% 3维：生成旋转矩阵,方向余弦矩阵x='x',y='y',z='z':check=1旋转，check=2方向余弦. check表达不正确
% 2020.11.12
r=zeros(3,3);
if(check==1)
    switch xyz
        case 'x'
            r=[1 0 0;
             0, cos(radian), -sin(radian);
             0, sin(radian), cos(radian)];
        case 'y'
          r=[cos(radian), 0, sin(radian);
           0,   1,0;
           -sin(radian),0 , cos(radian)];
         case 'z'
          r=[cos(radian), -sin(radian),0;
            sin(radian),cos(radian),0;
            0, 0, 1];
    end
elseif(check==2)
        switch xyz
         case 'x'
            r=[1 0 0;
             0, cos(radian), sin(radian);
             0, -sin(radian), cos(radian)];
        case 'y'
        r=[cos(radian), 0, -sin(radian);
           0,   1,0;
           sin(radian),0 , cos(radian)];
        case 'z'
        r=[cos(radian), sin(radian),0;
            -sin(radian),cos(radian),0;
            0, 0, 1];
    end
end
mrot=r;
end