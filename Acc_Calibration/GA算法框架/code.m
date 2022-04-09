function ret = code(lenchrome,LB,UB,flag)
%% Function Introduction
%线性编码产生初始种群
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Function realize  
    switch(flag)
        case 1
            ret = liner_code(lenchrome,LB,UB);%线性编码
        case 2
            ret = adaptive_code(lenchrome,LB,UB);%部分编码     
    end
end

function res = liner_code(lenchrome,LB,UB)
	flag=0;
    LB=LB';UB=UB';%将原本的列向量变为行向量
    while flag==0
		pick1 = rand(1,6);
		ret1=LB(1:6)+(UB(1:6)-LB(1:6)).*pick1;
		pick2 = 0.5*rand(1,3);
		ret2=LB(7:9)+(UB(7:9)-LB(7:9)).*pick2;
		ret = [ret1 ret2];
		flag=test(lenchrome,LB,UB,ret);%边界检测
    end
    res = ret;
end

function res = adaptive_code(lenchrome,LB,UB)
   flag=0;
    while flag==0
		pick = rand(1,lenchrome);
		ret  = LB' +(UB' - LB').*pick; 
		flag=test(lenchrome,LB,UB,ret);%边界检测
    end
    res = ret;
end

