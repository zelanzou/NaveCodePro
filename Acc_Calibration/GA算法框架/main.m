close all;                                           % close all figures
clc;                                                 % clear cmd text
clear ;                                              % clear all RAM
disp('加计标定');
popsize = 100;%种群大小
gensize = 500;%迭代次数
lenchrome = 9;%染色体个数即参数个数,9个
LB=[0.1,0.1,0.1,0.95,0.95,0.95,-0.06,-0.06,-0.06]';   % 参数下界
UB=[0.2,0.2,0.2,1.05,1.05,1.05,0.03,0.03,0.03]';   % 参数上界
individuals=struct('fitness',zeros(1,popsize),'chrome',[]);%种群结构体包含适应度和染色体，结构体更容易编写程序
bf=[];%缓存种群最佳适应度，记录每一代进化中朿好的适应度和平均适应度，可以用来绘图反映是否过早收敛
af=[];%缓存种群平均适应适应度
navs=[];

%产生初始种群
for i=1:popsize
    individuals.chrome(i,:) = code(lenchrome,LB,UB,1);%随机产生
    x=individuals.chrome(i,:);%随机种群的染色体（参数）
    individuals.fitness(i)=searchfun(x,0);%计算染色体适应度
end
%适应度变换
individuals.fitness = fitness_change(individuals,popsize,2);
%找最好的染色使
[bestfitness,bestindex]=max(individuals.fitness);%得到初始化种群的最大的适应度函数值和位置。
bestchrome=individuals.chrome(bestindex,:);%适应度最好的染色体
sumfitness = sum(individuals.fitness);
fmax = max(individuals.fitness);
favg = sumfitness/popsize; 
for i = 1:gensize
	num = i;
	individuals = select(popsize,individuals,sumfitness,1);%选择
	individuals.chrome = cross(popsize,individuals,lenchrome,LB,UB,num,5);%交叉
	individuals.chrome = mutation(popsize,individuals,fmax,favg,lenchrome,LB,UB,[1 gensize],num,5);%变异
    %计算适应度
    for j = 1:popsize
        x=individuals.chrome(j,:);
        individuals.fitness(j) = searchfun(x,0);
    end
	individuals.fitness = fitness_change(individuals,popsize,2);
    sumfitness = sum(individuals.fitness);
    fmax = max(individuals.fitness);
    favg = sumfitness/popsize; 
	[newbestfitness,newbestindex] = max(individuals.fitness);
	if bestfitness < newbestfitness%比较最优适应度
       bestfitness = newbestfitness; 
       bestchrome = individuals.chrome(newbestindex,:);%当前最优的参数  
	end
	navs=[navs;bestchrome];
    bf = [bf bestfitness];
    af = [af favg];
end
disp(bestchrome);%最后一个bestchrome
relative_erro = erro(bestchrome);
splot(bf,af,navs);