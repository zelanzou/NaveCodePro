function [bestchrome,x0] = GA_main(vsb)
%% Algorithm Introduction
%  A test model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author :zouzelan
% Version : V3.0 
% Date :
% File : 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

g0 = 9.794281725483827;
Sz = (vsb(1,3)-vsb(5,3))/(2*g0);
Sx = (vsb(9,1)-vsb(13,1))/(2*g0);
Sy = (vsb(17,2)-vsb(21,2))/(2*g0);
Bz = (vsb(1,3)+vsb(5,3))/2;
Bx = (vsb(9,1)+vsb(13,1))/2;
By = (vsb(17,2)+vsb(21,2))/2;

Ks=diag([Sx,Sy,Sz])^(-1);
x0 = [Bx,By,Bz,Ks(1,1),Ks(2,2),Ks(3,3),0,0,0];
Sx=Ks(1,1);  Sy=Ks(2,2);  Sz=Ks(3,3);
%% Algorithm time

%% Memory allocation
bf=[];%缓存种群最佳适应度，记录每一代进化中最好的适应度和平均适应度，可以用来绘图反映是否过早收敛
af=[];%缓存种群平均适应适应度
navs = [];
best = [];
fit = [];
chrome = [];
%% Initial data 
popsize = 50;%种群大小
gensize = 100;%迭代次数
lenchrome = 9;%染色体个数即参数个数,9个


% LB=[Bx-0.002,By-0.001,Bz-0.003,Sx-0.005,Sy-0.005,Sz-0.005,-0.02,-0.02,-0.02]';   % 参数下界
% UB=[Bx+0.001,By+0.001,Bz+0.003,Sx+0.005,Sy+0.005,Sz+0.005,0,0,0]';   % 参数上界

LB=[Bx-0.05,By-0.05,Bz-0.05,Sx-0.01,Sy-0.01,Sz-0.01,-0.02,-0.02,-0.02]';   % 参数下界  2020.09.03修改
UB=[Bx+0.05,By+0.05,Bz+0.05,Sx+0.01,Sy+0.01,Sz+0.01,0.02,0.02,0.02]';   % 参数上界




%% Function index
% buf=zeros(5*gensize,9);
% kk=1;
%% Algorithm develop
 for k = 1:5   %
     
	for i=1:popsize
         for m= 1:popsize
            individuals.chrome(m,:) = code(lenchrome,LB,UB,2);%随机产生
            x=individuals.chrome(m,:);%随机种群的染色体（参数）
            individuals.fitness(m)=searchfun(x,vsb,g0,2);%计算染色体适应度
         end
         [F,code_index] = sort(individuals.fitness);
         individuals.chrome(i,:) = individuals.chrome(code_index(1),:);
         individuals.fitness(i) = F(1);
    end
    sumfitness = sum(individuals.fitness);
	favg = sumfitness/popsize;
	fmin = min(individuals.fitness);
	alpha = favg/(favg-fmin);beta = fmin*favg/(favg-fmin);
	for i = 1:popsize
		fitness(i)= alpha*individuals.fitness(i)+beta;
    end
    individuals.fitness = fitness;
    
	%找最好的染色使
	[bestfitness,bestindex]=min(individuals.fitness);%得到初始化种群的最大的适应度函数值和位置。
	bestchrome=individuals.chrome(bestindex,:);%适应度最好的染色体
    individuals.fitness = 1./individuals.fitness;
	sumfitness = sum(individuals.fitness);
	fmax = max(individuals.fitness);
	favg = sumfitness/popsize; 

	for i = 1:gensize

		num = i;%计数
		individuals = select(popsize,individuals,sumfitness,1);%选择
		
		individuals.chrome = GA_cross(lenchrome,individuals,popsize,LB,UB,num,fmax,favg,5);%交叉
		
		individuals.chrome = GA_mutation(popsize,individuals,fmax,favg,lenchrome,LB,UB,[1 gensize],num,5);%变异
		
		%计算适应度
        for j = 1:popsize
            x=individuals.chrome(j,:);
            individuals.fitness(j)=searchfun(x,vsb,g0,2);%计算染色体适应度
        end
        
%         individuals = elitism_save(individuals);%精英保留策略
        
		[newbestfitness,newbestindex] = min(individuals.fitness);
        if bestfitness >=newbestfitness%比较最优适应度
           bestfitness = newbestfitness; 
           bestchrome = individuals.chrome(newbestindex,:);%当前最优的参数  
        end
      
        sumfitness = sum(individuals.fitness);
        favg = sumfitness/popsize;
		bf = [bf bestfitness];
		af = [af favg];
		navs=[navs;bestchrome];

        
        individuals.fitness = 1./individuals.fitness;
        sumfitness = sum(individuals.fitness);
        fmax = max(individuals.fitness);
        favg = sumfitness/popsize; 
        
%         buf(kk,:)=bestchrome;  kk=kk+1;
    end
    [fitness,index] = min(bf);
    chrome = navs(index(1),:);
%     % mei10
%     plot(buf(1:kk-1,4));
%      Dd=Dd*0.6;
%      LB=(chrome-Dd)';
%      UB=(chrome+Dd)';
 end
bestchrome = chrome;

end


