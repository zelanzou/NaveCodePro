% Note : This main is data process,and kf parameter setting, Users could
% use gps ,imu, imuerr,dap0 and avp in other function.
% Example :res= kf(imu,gps,davp0,imuerr,avp0)
clc;
clear;
%% data1
load('trjNm_R11.mat')
imu_adis=[trjN_R1.adis(:,1:3)*pi/180,trjN_R1.adis(:,4:6),trjN_R1.tad'];
gnss_ref = [trjN_R1.ngnss(1:2:end,1:2)*pi/180,trjN_R1.ngnss(1:2:end,3),repmat(1,length(trjN_R1.ngnss(1:2:end,1:2)),1),trjN_R1.ngnss(1:2:end,4)]; %2HZ
avpm = interp1(trjN_R1.tref, trjN_R1.avpm(:,1:9),trjN_R1.tad,'linear');  %线性插值
trj.avp = [avpm(:,1:3)*pi/180,avpm(:,4:6),avpm(:,7:8)*pi/180,avpm(:,9),trjN_R1.tad'];
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
eb =[143.7;100.2;145.4]*glv.dph;
db = [3000;3000;3000]*glv.ug;
web = 1*[0.23;0.18;0.23]*glv.dpsh;
wdb = 1*[56;52;66]*glv.ugpsHz; 
fb0 = mean(imu_adis(1:8*100,4:6));
wb0 = mean(imu_adis(1:8*100,1:3));

gps = [zeros(length(gnss_ref),3),gnss_ref];
imu = [imu_adis(:,1:3) - wb0,imu_adis(:,4:6),imu_adis(:,end)];
imuerr = struct('eb',eb,'db',db,'web',web,'wdb',wdb);
davp0 = avpseterr([1800;1800;60], [0.05;0.05;0.05], [3;3;3]); %初始误差 0.5°
avp0 = [trj.avp(1,1:6)';gnss_ref(1,1:3)'];

%% data2
load('2021年05月31日18时45分12秒_KVH_DATA.txt');
load('2021年05月31日18时45分12秒_SPAN_DATA.txt');
A =  X2021_05_31_18_45_12__SPAN_DATA; 
B =  X2021_05_31_18_45_12__KVH_DATA;
trj.avp = [A(:,11),A(:,10),A(:,12),A(:,8),A(:,7),A(:,9),A(:,4:6),A(:,3)];
trj.avp = [trj.avp(:,1:3)*pi/180, trj.avp(:,4:6), trj.avp(:,7:8).*pi/180, trj.avp(:,9:10)];
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
imuerr = imuerrset(10, 1000, 0.2, 100);
davp = avpseterr([-30;30;30], [0.1;0.1;0.1]*1, [1;1;1]*2);              % 设定初始姿态、速度、位置误差

mavp  = trj.avp;
for i=1:length(mavp)
    if(mavp(i,3)>pi )
        mavp(i,3) = 2*pi - mavp(i,3);
    else
        mavp(i,3) = - mavp(i,3);
    end
end
imu = [B(:,4:6)*pi/180,B(:,7:9),B(:,3)];
avp0 = mavp(1,:)';
%res  = SINS186Att_Vel_TransferAlignment(imu, mavp, imuerr, davp,avp0);
%% data3
load('0604SPAN.mat');
load('0604KVH.mat');
load('0604DSP.mat');
trj.avp = [SPAN(:,11),SPAN(:,10),SPAN(:,12),SPAN(:,8),SPAN(:,7),SPAN(:,9),SPAN(:,4:6),SPAN(:,3)];
trj.avp = [trj.avp(:,1:3)*pi/180, trj.avp(:,4:6), trj.avp(:,7:8).*pi/180, trj.avp(:,9:10)];
global glv 
glv = globalParameter(trj.avp(1,7),trj.avp(1,9));
imuerr = imuerrset(10, 1000, 0.2, 100);
davp = avpseterr([-60;60;300], [1;1;1]*1, [5;5;5]*2);              % 设定初始姿态、速度、位置误差

mavp  = trj.avp;
for i=1:length(mavp)
    if(mavp(i,3)>pi )
        mavp(i,3) = 2*pi - mavp(i,3);
    else
        mavp(i,3) = - mavp(i,3);
    end
end
mavp(:,end) = mavp(:,end) -mavp(1,end) ;
avp0 = mavp(1,:)';

imukvh = [KVH(:,4:6)*pi/180,KVH(:,7:9),KVH(:,3)-KVH(1,3)];
imustim = [ DSP(:,2:4)*glv.deg DSP(:,5:7) imukvh(:,end)];
%res  = SINS186Att_Vel_TransferAlignment(imukvh, mavp, imuerr, davp,avp0);