function imu  = WaveDenoise(imu) 
s1 = imu(:,1);s2 = imu(:,2);s3 = imu(:,3);
s4 = imu(:,4);s5 = imu(:,5);s6 = imu(:,6);
% thrrr3='heursure';
wavec='sym4';
M3 = 5 ;%
xds1=wden(s1,'heursure','s','mln',M3,wavec);
xds2=wden(s2,'heursure','s','mln',M3,wavec);
xds3=wden(s3,'heursure','s','mln',M3,wavec);
xds4=wden(s4,'heursure','s','mln',M3,wavec);
xds5=wden(s5,'heursure','s','mln',M3,wavec);
xds6=wden(s6,'heursure','s','mln',M3,wavec);
imu = [xds1 xds2 xds3 xds4 xds5 xds6 imu(:,7)];
end