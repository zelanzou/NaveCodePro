function [imu, avp0] = Myavp2imu(avp)
% Simulate SIMU sensor outputs from attitude, velocity & position profile.
%
% Prototype: [imu, avp0] = avp2imu(avp)
% Input: avp = [att,vn,pos,t]
% Outputs: imu = [wm,vm,t]  增量形式
%          avp0 = init [att,vn,pos]
%
% See also  ap2avp, ap2imu, trajsimu, insupdate.

% Copyright(c) 2009-2014, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 15/10/2013, 15/03/2014
global glv
    if size(avp,2)<9, avp = ap2avp(avp, diff(avp(1:2,end))); end
    len = size(avp,1);  ts = avp(2,10)-avp(1,10);  ts2 = ts/2;
    Cbn_1 = a2mat(avp(1,1:3)')';  vn_1 = avp(1,4:6)';  pos_1 = avp(1,7:9)';
    wm_1 = [0;0;0];  vm_1 = [0;0;0];
    imu = zeros(len, 6);
    timebar(1, len, 'Trajectory inversion avp->imu.');
    for k=2:len  % begin from 2
        Cnb = a2mat(avp(k,1:3));   vn = avp(k,4:6)';   pos = avp(k,7:9)';
        eth =  EarthParameter((pos_1+pos)/2, (vn_1+vn)/2);
%         wbts = m2rv(Cbn_1*Cnb);
%         phim = wbts + (Cbn_1+Cnb')*(eth.wnin*ts2);
        phim = m2rv(Cbn_1*rv2m(eth.wnin*ts)*Cnb);
        wm = (glv.I33+askew(wm_1/12))^-1*phim; % using previous subsample: phim = wm + 1/12*cross(wm_1,wm)
        dvbm = Cbn_1*qmulv(rv2q(eth.wnin*ts2), vn-vn_1-eth.gcc*ts); % sins
%         vm = (glv.I33+askew(wm/2))^-1*dvbm;  % dvbm = vm + 1/2*cross(wm,vm)
        vm = (glv.I33+askew(wm/2+wm_1/12))^-1*(dvbm-cros(vm_1,wm)/12);  % dvbm = vm + 1/2*cross(wm,vm)
        imu(k,:) = [wm; vm]';
        Cbn_1 = Cnb'; vn_1 = vn; pos_1 = pos; wm_1 = wm; vm_1 = vm;
        timebar;
    end
    imu = [imu(2:end,:), avp(2:end,10)];
    avp0 = avp(1,1:9)';

