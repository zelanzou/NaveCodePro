function module_value = static_pos(x,fs,flag)
% g0 = 9.794731130106427;
    Ka = [x(4),  0,   0
          x(7), x(5), 0
          x(8), x(9), x(6)];
    bias = [x(1) x(2) x(3)];
    fs = fs';
     switch(flag)
		case 1 %模型1
			fb = Ka*fs;
			fb(1,:)=fb(1,:)-bias(1); fb(2,:)=fb(2,:)-bias(2); fb(3,:)=fb(3,:)-bias(3);
		case 2 %模型2
			fb = Ka* [fs(1,:)-bias(1);fs(2,:)-bias(2);fs(3,:)-bias(3)];
        case 3 %模型3
            fb = fs+Ka*fs;
            fb(1,:)=fb(1,:)-bias(1); fb(2,:)=fb(2,:)-bias(2); fb(3,:)=fb(3,:)-bias(3);
     end
    fb = fb';
    %计算模值
     for j = 1:length(fb)
        value(j) = sqrt(sum(fb(j,:).^2));
     end
     module_value = value;
end
