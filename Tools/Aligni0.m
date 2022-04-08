function [att0, res] = Aligni0(imu, pos, ts)
%严恭敏的利用gps初始对准程序
global glv
    if nargin<3,  ts = imu(2,end)-imu(1,end);  end
    nn = 2; nts = nn*ts;  ratio = 1; % 0.995;
    len = fix(length(imu)/nn)*nn;
    if length(pos)>4 , pos=pos(4:6); 
    elseif length(pos)==1, pos=[pos;0;0]; end
    eth = EarthParameter(pos);  lat = pos(1);  g0 = -eth.gn(3);
    qib0b = [1; 0; 0; 0];
    [vib0, vi0, pib0, pi0, vib0_1, vi0_1] = setvals(zeros(3,1));
    [pib0k, pi0k, vi0k, vib0k, fi0k, fib0k, attk, attkv] = prealloc(len/nn, 3);
    k0 = fix(5/ts); % exculde the first 5s
    ki = timebar(nn, len, 'Initial align based on inertial frame.');
    for k=1:nn:len-nn+1
        wvm = imu(k:k+nn-1, 1:6);  kts = (k+nn-1)*ts;
        [phim, dvbm] = cnscl(wvm);
        fib0 = qmulv(qib0b, dvbm)/nts;   % f
        vib0 = vib0 + fib0*nts;          % vel
        pib0 = ratio*pib0 + (vib0_1+vib0)*nts/2;  vib0_1 = vib0; % pos
%         fi0 = [eth.cl*cos(kts*glv.wie);eth.cl*sin(kts*glv.wie);eth.sl]*g0;
%         vi0 = vi0 + fi0*nts;
%         pi0 = ratio*pi0 + (vi0_1+vi0)*nts/2;      vi0_1 = vi0;
        [fi0, vi0, pi0] = i0fvp(kts, lat);
        qib0b = qupdt(qib0b, phim);  % qib0b updating
        pib0k(ki,:) = pib0'; vib0k(ki,:) = vib0'; fib0k(ki,:) = fib0'; % recording
        pi0k(ki,:) = pi0';   vi0k(ki,:) = vi0';   fi0k(ki,:) = fi0';
        if k>k0
            k1 = fix(ki/2);
            swiet = sin(kts*glv.wie); cwiet = cos(kts*glv.wie);
            Cni0 = [-swiet,cwiet,0; 
                -eth.sl*cwiet,-eth.sl*swiet,eth.cl; 
                eth.cl*cwiet,eth.cl*swiet,eth.sl];
            qni0 = m2qua(Cni0);
            qi0ib0 = dv2atti(vi0k(k1,:)', vi0, vib0k(k1,:)', vib0);
            qnb = qmul(qmul(qni0,qi0ib0),qib0b);
            attkv(ki,:) = q2att(qnb)';    % using vel
            qi0ib0 = dv2atti(pi0k(k1,:)', pi0, pib0k(k1,:)', pib0);
            qnb = qmul(qmul(qni0,qi0ib0),qib0b);
            attk(ki,:) = q2att(qnb)';     % using pos
       end
       ki = timebar;
    end
    k0 = fix(k0/nn)+1;
%     attk(1:k0,:) = repmat(attk(k0+1,:),k0,1);
    Cni0 = [0,1,0; -eth.sl,0,eth.cl;  eth.cl,0,eth.sl];
    att0 = q2att(qmul(m2qua(Cni0),qi0ib0));
    attk(1:k0,:) = repmat(att0',k0,1);
    attkv(1:k0,:) = repmat(attkv(k0+1,:),k0,1);
    tk = imu(nn:nn:length(attk)*nn,7); attk(:,4) = tk; attkv(:,4) = tk;
    res = varpack(lat, nts, vib0k, pib0k, fib0k, vi0k, pi0k, fi0k, attk, attkv, att0); 
    att0 = attk(end,1:3)';
    resdisp('Initial align attitudes (arcdeg)', att0/glv.deg);
    ai0plot(nts, attk, attkv);
    
function ai0plot(ts, attk, attkv)
global glv
    t = (1:length(attk))'*ts;
    myfigure;
    subplot(211), plot(t, attk(:,1:2)/glv.deg), xygo('pr');
        hold on,  plot(t, attkv(:,1:2)/glv.deg, 'm:'),
    subplot(212), plot(t, attk(:,3)/glv.deg), xygo('y');
        hold on,  plot(t, attkv(:,3)/glv.deg, 'm:'), legend('i0 pos', 'i0 vel'); title(sprintf('\\psi=%.4f', attk(end,3)/glv.deg));
