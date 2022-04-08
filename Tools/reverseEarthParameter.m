function eth = reverseEarthParameter(pos,vn)
Re = 6378137;
f = 1/298.257;
wie = -7.2921151467e-5;
g0 = 9.7803267714; %赤道附近的重力
eth = [];
eth.pos = pos;  eth.vn = vn;
eth.sl = sin(pos(1));  eth.cl = cos(pos(1));  eth.tl = eth.sl/eth.cl; %pos(1)就是L
eth.sl2 = eth.sl*eth.sl;  sl4 = eth.sl2*eth.sl2;
eth.RMh = Re*(1-2*f+3*f*eth.sl2);%zzl 2020.10.20
eth.RNh = Re*(1-f*eth.sl2); eth.clRNh = eth.cl*eth.RNh;
eth.wnie = [0; wie*eth.cl; wie*eth.sl];
vE_RNh = vn(1)/eth.RNh;
eth.wnen = [-vn(2)/eth.RMh; vE_RNh; vE_RNh*eth.tl];
eth.wnin = eth.wnie + eth.wnen;
eth.wnien = eth.wnie + eth.wnin;
eth.g = g0*(1+5.27094e-3*eth.sl2+2.32718e-5*sl4)-3.086e-6*pos(3); %当地重力
eth.gn = [0;0;-eth.g];
eth.gcc =  [ eth.wnien(3)*vn(2)-eth.wnien(2)*vn(3);  % faster than previous line 有害加速度
             eth.wnien(1)*vn(3)-eth.wnien(3)*vn(1);
             eth.wnien(2)*vn(1)-eth.wnien(1)*vn(2)+eth.gn(3) ];
end