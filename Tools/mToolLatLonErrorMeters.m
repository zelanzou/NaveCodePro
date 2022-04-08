function err=mToolLatLonErrorMeters(rflat,rflon,lat,lon,rfH)
% 工具：转换经纬度误差 为米，行向量
%%======================================================================= 
%input:
%   lat:解算的纬度 1*N
%   lon:解算的经度 1*N
%   relat:参考轨迹的纬度 1*N
%   reflon:参考轨迹的经度 1*N
%   rfH: 参考高度
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author : xutongxu
% Version :  
% Date : 2020.9.8
% File : 
%output :err.erx:纬度误差
%             err.ery:经度误差
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

N=length(rflat);
erlat=lat-rflat;
erlon=lon-rflon;
Re=6378137.0;
f=0.003352813177897;

Rm=zeros(1,N);   Rn=Rm;
ery=zeros(1,N);   erx=ery;
if(nargin>4)
	for i=1:N
		L0=rflat(i);
		Rm(i)=Re*(1-2*f+3*f*sin(L0)^2) +rfH(i);
		Rn(i)=Re*(1+f*sin(L0)^2) +rfH(i); 
		ery(i)=Rm(i)*erlat(i);
		erx(i)=Rn(i).*erlon(i)*cos(L0);
	end
else
    for i=1:N
		L0=rflat(i);
		Rm(i)=Re*(1-2*f+3*f*sin(L0)^2);
		Rn(i)=Re*(1+f*sin(L0)^2); 
		ery(i)=Rm(i)*erlat(i);
		erx(i)=Rn(i).*erlon(i)*cos(L0);
    end
end

% theta*R=arc
% erN=Rm.*erlat;
% erE=Rn.*erlon;
err=varpack(ery,erx);
end