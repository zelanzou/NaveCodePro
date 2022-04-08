function [att,vel,pos,qua,Cnb,eth] = avp_update(wbib,fb,Cnb,qua,p0,v0,dt)
%纯惯导解算
eth = EarthParameter(p0,v0);
fn = Cnb*fb;
an = fn + eth.gcc;
vel = v0 + an*dt;

Mpv = [0 1/eth.RMh 0;1/eth.clRNh 0 0;0 0 1];
pos = p0 + Mpv*vel*dt;

eth = EarthParameter(pos,vel); %这里是否更新地球参数？
wbnb = wbib - Cnb'*eth.wnin;%b系下
qua = qupdt(qua,wbnb*dt);%wbnb*dt是相当于等效旋转增量rv
Cnb = q2mat(qua);
att = m2att(Cnb);
end