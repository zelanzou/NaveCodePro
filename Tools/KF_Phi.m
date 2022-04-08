%% This is a 15-dimensional system state discrete model

%Note£º 1, the input parameter "flag" is dimension
%            2, the eth should use the tools function of EarthParameter
%---------------------------------------------------------------------------%
function F = KF_Phi(eth,Cnb,fb,flag)
%% Ò»²½×´Ì¬×ªÒÆ¾ØÕó
Maa = -[0,            -eth.wnin(3), eth.wnin(2);
       eth.wnin(3),  0,            -eth.wnin(1);
       -eth.wnin(2), eth.wnin(1),  0           ];% Maa = -askew(ins.eth.wnin);
   
Mav = [0,              -1/eth.RMh, 0
       1/eth.RNh,      0,          0;
       eth.tl/eth.RNh, 0,          0];
M1 = [ 0,            0, 0;
       -eth.wnie(3), 0, 0;
        eth.wnie(2), 0, 0 ];
M2 = [0,                            0, eth.vn(2)/(eth.RMh)^2;
      0,                            0, -eth.vn(1)/(eth.RNh)^2;
      eth.vn(1)/(eth.clRNh*eth.cl), 0, -eth.vn(1)*eth.tl/(eth.RNh)^2];
Map = M1+M2;
fn = Cnb *fb;
Mva = askew(fn);
Mvv = askew(eth.vn)*Mav - askew(eth.wnien);
Mvp = askew(eth.vn)*(2*M1+M2);
scl = eth.sl*eth.cl; g0 = 9.7803267714;
Mvp(3,1) = Mvp(3,1)-g0*(5.27094e-3*2*scl+2.32718e-5*4*eth.sl2*scl);
Mvp(3,3) = Mvp(3,3)+3.086e-6;

Mpv = [0,                 1/eth.RMh, 0;
       1/(eth.RNh*eth.cl),0,         0;
       0,                 0,         1];
Mpp = [0,                                0, -eth.vn(2)/(eth.RMh)^2;
       eth.vn(1)*eth.tl/(eth.RNh*eth.cl),0, -eth.vn(1)/((eth.RNh)^2*eth.cl);
       0,                                0, 0                             ];
O33 = zeros(3);
	%%   euler dvn   dpos   eb     db
	F = [ Maa    Mav    Map    -Cnb     O33 
          Mva    Mvv    Mvp     O33     Cnb 
          O33    Mpv    Mpp     O33     O33
          zeros(6,15) ];
    if (flag ~= 15)
        F(flag,flag) = 0;
    end
        
end