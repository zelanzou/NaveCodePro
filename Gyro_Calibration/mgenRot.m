function mrot=mgenRot(radian,xyz,check)
% 3维：生成旋转矩阵,方向余弦矩阵x='x',y='y',z='z':check=1旋转，check=2方向余弦.
% 函数原型：mrot=mgenRot(radian,xyz,check)
r=zeros(3,3);

if(check==1)
    switch xyz
        case 'x'
            r=[1 0 0;
             0, cos(radian), -sin(radian);
             0, sin(radian), cos(radian)];
        case 'y'
          r=[cos(radian), 0, sin(radian);
           0,   1,0;
           -sin(radian),0 , cos(radian)];
         case 'z'
          r=[cos(radian), -sin(radian),0;
            sin(radian),cos(radian),0;
            0, 0, 1];
    end
elseif(check==2)
        switch xyz
         case 'x'
            r=[1 0 0;
             0, cos(radian), sin(radian);
             0, -sin(radian), cos(radian)];
        case 'y'
        r=[cos(radian), 0, -sin(radian);
           0,   1,0;
           sin(radian),0 , cos(radian)];
        case 'z'
        r=[cos(radian), sin(radian),0;
            -sin(radian),cos(radian),0;
            0, 0, 1];
    end
end
mrot=r;
return 