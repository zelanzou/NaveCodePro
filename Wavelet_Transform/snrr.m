function z=snrr(x,y)%计算信噪比，x是原始信号，y是去噪后信号
y1=sum(x.^2);
y2=sum((y-x).^2);
z=10*log10((y1/y2));
end
