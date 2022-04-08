function res=Arw_Vrw2std(arw,vrw,dt,pos0,flag)
% 严恭敏程序中 ARW, VRW 与 噪声标准差的转换
% 函数原型： res=Arw_Vrw2std(arw,vrw,dt,pos0,flag)
% 输入参数：
%    dt：采样间隔，单位 s.    pos0:[纬度; 经度; 高度], 纬度、经度单位：rad，高度单位m
%  falg=1时, arw单位 度/sqrt(h),  vrw单位 ug/sqrt(Hz)
%  flag=2时, arw单位 度/s,  vrw单位 m/s2, 为静止数据标准差.
%
%  输出参数;res: {arw, vrw, acc_std, gyro_std}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : Xu Tongxu
% Date : 2021.3.31
% File : 2015.NaveGo: a simulation framework for low-cost integrated navigation systems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% eth=earth(pos0,zeros(3,1));
g0  = mxGetGravity(pos0(1),pos0(3));
ug=g0*10^(-6);

if(flag==1)
    % 度/sqrt(h) 转 度/s
    gyro_std=(arw/(60*sqrt(dt)));         % 采样间隔 dt下，静态数据标准差，噪声
    % ug/sqrt(Hz) 转 m/s2
    acc_std=(vrw/sqrt(dt))*ug;            % 采样间隔 dt下，静态数据标准差，噪声
else
    if(flag==2)
        gyro_std=arw;     acc_std=vrw;
         %  度/s  转  度/sqrt(h)
         arw=gyro_std*60*sqrt(dt);
         
         % m/s2   转  ug/sqrt(Hz)
         vrw=acc_std*sqrt(dt)/ug;
    end
end

res=varpack(arw,vrw,gyro_std,acc_std);

end