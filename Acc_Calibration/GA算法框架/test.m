function flag=test(lenchrome,LB,UB,code)
%检测越界函数
%lenchrom input: 染色体长度,也就是一条染色体上的基因个数
%bound input :变量的取值范围
%code input :染色体的编码值
flag=1;
[m,n]=size(code);
for i=1:n         
    if code(i)<LB(i)||code(i)>UB(i)%越界，重新编码
        flag=0;
    end
end