function kf = ukf_filter(kf,avp,imu,gamma,Wm,Wc,nts,flag,Z,H)
    m = size(kf.Pk,1);
%     q = size(kf.Qk,1);
    % r = size(R,1);
%     n = m+q;
    n = m;
    switch m
        case 15
            if (flag == 1)
%                 Paug = [kf.Pk,zeros(m,q);zeros(q,m),kf.Qk];
                Paug = [kf.Pk];
                L = chol( Paug,'lower');
%                 Xaug = [kf.Xk;zeros(q,1)];
                Xaug = [kf.Xk];
                x = repmat(Xaug,1,n);
                xsigma0 = [Xaug,x+gamma*L,x-gamma*L];%15*31
                for i = 1:(2*n+1)
%                     Xsigma1(:,i) = state_function(xsigma0(1:15,i),xsigma0(16:21,i),imu,avp,nts); %15*31
                       Xsigma1(:,i) = state_function(xsigma0(1:15,i),imu,avp,nts); %15*31
                end
                %时间更新
                X0 =Xsigma1* Wm ; %15*1
                for i = 1:(2*n+1)
                    xerr(:,i) = Xsigma1(:,i)-X0; %15*1  
                end
                
                kf.Qk = kf.Qt*nts;
                P0= xerr*Wc*xerr'+ kf.Gammak*kf.Qk*kf.Gammak'; 
                kf.Xkk_1 = X0;
                kf.Pkk_1 = P0;
                %量测更新
            %	[Z,H] = measure_function();
                kf.Hk = H;
                kf.PXZkk_1 = kf.Pkk_1*kf.Hk';%状态一步预测与量测一步预测的协均方误差阵
                kf.rk = Z - kf.Hk*kf.Xkk_1;%残差
                kf.PZZkk_1 = kf.Hk*kf.PXZkk_1 + kf.Rk;%量测一步预测均方误差阵
                kf.Kk = kf.PXZkk_1*invbc(kf.PZZkk_1);
                kf.Xk = kf.Xkk_1 + kf.Kk*kf.rk;
                kf.Pk = kf.Pkk_1 - kf.Kk*kf.PZZkk_1*kf.Kk';
                kf.Pk = (kf.Pk+kf.Pk')/2;
            elseif (flag == 0)
%                 Paug = [kf.Pk,zeros(m,q);zeros(q,m),kf.Qk];
                Paug = [kf.Pk];
                L = chol( Paug,'lower');
%                 Xaug = [kf.Xk;zeros(q,1)];
                Xaug = [kf.Xk];
                x = repmat(Xaug,1,n);
                xsigma0 = [Xaug,x+gamma*L,x-gamma*L];%15*31
                for i = 1:(2*n+1)
%                     Xsigma1(:,i) = state_function(xsigma0(1:15,i),xsigma0(16:21,i),imu,avp,nts); %15*31
                      Xsigma1(:,i) = state_function(xsigma0(1:15,i),imu,avp,nts); %15*31
                end
                %时间更新
                X0 =Xsigma1* Wm ; %15*1
                for i = 1:(2*n+1)
                    xerr(:,i) = Xsigma1(:,i)-X0; %15*1  
                end
                kf.Qk = kf.Qt*nts;
                P0= xerr*Wc*xerr' + kf.Gammak*kf.Qk*kf.Gammak';
                kf.Xk = X0;
                kf.Pk = P0;
            end
        otherwise
    end
end