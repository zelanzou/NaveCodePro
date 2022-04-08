function gravity=mxGetGravity(L,h)
% 根据纬度、高度 计算重力加速度
% 函数原型：gravity=mxGetGravity(L,h)
% L:纬度，单位：弧度 h:高度，相对海平面
% 参考文献：《惯性导航》第二版，秦永元，178页，公式（7.2.17）,7.2.18
Re=6378137;% m
var1=(1+0.00193185138639*sin(L)^2);
var2=sqrt(1-0.00669437999013*sin(L)^2);
g=978.03267714*var1/var2;
gra=g*Re^2/((Re+h)^2); % cm/s^2
gravity=gra/100;
return ;