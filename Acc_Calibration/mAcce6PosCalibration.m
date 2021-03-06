function res=mAcce6PosCalibration(vxu1,vxd1,vyu4,vyd2,vzu1,vzd3,g0)
% 6 位置修正标定法
%   vector_s=Sa*Asba*vector_b + bisv; 
%  函数原型： res=mAcce6PosCalibration(vxu1,vxd1,vyu4,vyd2,vzu1,vzd3,g0)
%---------------------------------------------------------------------------------------------------------%
% 输入变量：
%        静态12位位置平均值。单位：rad/s
% 输出变量：
%        res:结构体
%        { SFm-修正,    SF-未修正(up-down),  bisv-(零偏电压), sbis -(未修正零偏),  Asbm- 非正交矩阵}
%---------------------------------------------------------------------------------------------------------%
% Author : xutongxu
% Version : 1
% Date : 2020.2.18
% File : Xu T, Xu X, Xu D, et al. A Novel Calibration Method Using Six Positions for MEMS Triaxial Accelerometer[J]. IEEE Transactions on Instrumentation and Measurement, 2020, 70: 1-11. 
%---------------------------------------------------------------------------------------------------------%
 
Xt=zeros(10,12);
X0=zeros(1,15);
%%  fprintf('倾斜角beta : %.5f\n',betas*180/pi);
% 开始迭代计算
k=1; 
beta=0; qxy=0; qxz=0; qyx=0; qyz=0; qzx=0; qzy=0;
while(k<10) % 修正计算
    sfx=(vxu1(1)-vxd1(1))/(2*g0*cos(beta)*cos(qxy)*cos(qxz));  
    bisx=(vxu1(1)+vxd1(1))/2-g0*sfx*sin(beta)*sin(qxz);
    
    sfy=(vyu4(2)-vyd2(2))/(2*g0*cos(beta)*cos(qyx)*cos(qyz));  
    bisy=(vyu4(2)+vyd2(2))/2-g0*sfy*sin(beta)*sin(qyz)*cos(qyx);
    
    sfz=(vzu1(3)-vzd3(3))/(2*g0*cos(beta)*cos(qzx)*cos(qzy)); 
    bisz=(vzu1(3)+vzd3(3))/2+g0*sfz*sin(beta)*sin(qzx)*cos(qzy);
    
    qyx=asin((vzu1(2)-vzd3(2))/(2*g0*sfy*cos(beta)));
    qxz=asin((vyu4(1)-vyd2(1))/(2*g0*sfx*cos(beta)));
    qzy=asin((vxu1(3)-vxd1(3))/(2*g0*sfz*cos(beta)));

    qyz=asin((vxu1(2)-vxd1(2))/(-2*g0*cos(qyx)*sfy*cos(beta)));
    qxy=asin((vzu1(1)-vzd3(1))/(-2*g0*cos(qxz)*sfx*cos(beta)));
    qzx=asin((vyu4(3)-vyd2(3))/(-2*g0*cos(qzy)*sfz*cos(beta)));

    sbb1=(vxu1(2)+vxd1(2)-2*bisy)/( 2*g0*sfy*cos(qyx)*cos(qyz));
    sbb2=(vyu4(1)+vyd2(1)-2*bisx)/(-2*g0*sfx*cos(qxy)*cos(qxz));
    sbb3=(vzu1(2)+vzd3(2)-2*bisy)/( 2*g0*sfy*cos(qyx)*cos(qyz));

    beta=asin((sbb1+sbb2+sbb3)/3);
    X1=[sfx,sfy,sfz,bisx,bisy,bisz,qxy,qxz,qyx,qyz,qzx,qzy];
    Xt(k,:)=X1;
    if(k>1)
        err=abs(Xt(k,:)-Xt(k-1,:));
        if(max(err(1:12))<10^(-6))
            break;
        end
    else
        X0=X1;
    end
    k=k+1;
end
Xi=Xt(1:k,:);
%----修正标度因数、零偏
sfxm=X1(1);
sfym=X1(2);
sfzm=X1(3);
bisxm=X1(4);
bisym=X1(5);
biszm=X1(6);

SF=diag(X0(1:3));  SFm=diag([sfxm,sfym,sfzm]);
sbis=X0(4:6)';     bisv=[bisxm;bisym;biszm];  
%% vzu1 as 参考 :2020.2.9
% bisv=[vzu1(1)-sfxm*g0*(sin(beta)*sin(qxz)-cos(beta)*cos(qxz)*sin(qxy))
%       vzu1(2)-sfym*g0*(cos(beta)*sin(qyx)+sin(beta)*cos(qyx)*cos(qyz))
%       biszm];
% a=(vzu1(1)-bisv(1))/(sfxm*g0)-sin(beta)*sin(qxz);  qxy=asin(-a/(cos(beta)*cos(qxz)));

%%
Csx_b=mgenRot(qxz,'z',2)*mgenRot(qxy,'y',2);
Rx=Csx_b(1,:);                               % b系矢量在 s 系 x轴上投影
Csy_b=mgenRot(qyx,'x',2)*mgenRot(qyz,'z',2);
Ry=Csy_b(2,:);                               % b系矢量在 s 系 x轴上投影
Csz_b=mgenRot(qzy,'y',2)*mgenRot(qzx,'x',2);
Rz=Csz_b(3,:);                               % b系矢量在 s 系 x轴上投影
Asbm=[Rx;Ry;Rz];                              % 正交系b系矢量，在s系投影 变换矩阵
%% 添加 于：2020.5.30
ex=Csx_b*[1;0;0];  ey=Csy_b*[0;1;0];  ez=Csz_b*[0;0;1];
Ez=cross(ex,ey);  Ey=cross(Ez,ex);
Ex=ex/norm(ex);  Ey=Ey/norm(Ey);  Ez=Ez/norm(Ez);
CbbL=[Ex,Ey,Ez];
AsbL=Asbm*(CbbL'); % 12参数，转换到下三角模型
%%
res=varpack(SFm,SF,bisv,sbis,Asbm,Xi,k,beta,AsbL);

end