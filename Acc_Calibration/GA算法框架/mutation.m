function ret = mutation(popsize,individuals,fmax,favg,lenchrome,LB,UB,pop,num,flag3)
	for i = 1:popsize
		pick=rand;%随机产生变异概率
		while pick==0
			pick=rand;
		end
		index=ceil(pick*popsize);%取整被用来变异的染色体的编号
		switch(flag3)
			case 1 
				pm = SGA_mutation();%简单变异算子
			case 2 
				pm =  AGA_mutation(individuals,fmax,favg,index);%自适应变异算子
			case 3 
				pm = IAGA_mutation(individuals,fmax,favg,index);%线性自适应变异算子
			case 4 
				pm = HIAGA_mutation(individuals,fmax,favg,index);%异性自适应变异算子
			case 5 
				pm = TIAGA_mutation(individuals,fmax,favg,index,num);%迭代自适应变异算子	
		end
        pick=0.1*rand;%产生一个随机概率 pm在0.01-0.05之间，直接用rand太大，导致种群基本无变异
		while pick==0
			pick=rand;
		end
		if pick >pm      %如果产生的随机概率大于变异概率，那么此轮循环不变异，类似于交叉，变异是小概率事件
			continue;
		end
		flag=0;
		while flag==0
			pick=rand;%随机产生不为0的随机变异位置
			while pick==0
				pick=rand;
			end
			pos=ceil(pick*lenchrome);%产生变异的位置，也就是染色体上的第几个基因变异
			v=individuals.chrome(i,pos);%第i条染色体中pos基因
			v1=v-LB(pos);
			v2=UB(pos)-v;
			pick=rand;%开始变异，实值变异法
			if pick>0.5
				delta=v2*(1-pick^((1-pop(1)/pop(2))^2));
				individuals.chrome(i,pos)=v+delta;
			else
				delta=v1*(1-pick^((1-pop(1)/pop(2))^2));
				individuals.chrome(i,pos)=v-delta;
			end %变异结束			
			flag=test(lenchrome,LB,UB,individuals.chrome(i,:));
		end
	end
	ret=individuals.chrome;
end

function res  = SGA_mutation()
	pm = 0.1;%变异概率
	res = pm;
end

function res  = AGA_mutation(individuals,fmax,favg,index)
	k3=0.1;k4=0.1; 
	f = individuals.fitness(index);
	if f >= favg
        pm = k3*(fmax-f)/(fmax-favg);
    else
        pm = k4;
    end
	res = pm;
end

function res  = IAGA_mutation(individuals,fmax,favg,index)
	pm1=0.1; pm2=0.01;
	f = individuals.fitness(index);
	if f >= favg
        pm = pm1-(pm1-pm2)*(fmax-f)/(fmax-favg);
    else
        pm = pm1;
    end
	res = pm;
end

function res  = HIAGA_mutation(individuals,fmax,favg,index)
	pmmax = 0.1;pmmin = 0.006;beta = 5.512;
	f = individuals.fitness(index);
	alpha = (f - favg)/(fmax-favg);
    if f >= favg
        pm = 4*(pmmax-pmmin)/(exp(2*beta*alpha)+exp(-2*beta*alpha)+2)+pmmin;
    else
        pm = pmmax;
    end
	res = pm;
end

function res  = TIAGA_mutation(individuals,fmax,favg,index,num)
	phi = 0.1; pm2 = 0.01;
	f = individuals.fitness(index);
	alpha = (f - favg)/(fmax-favg);
	pm1 = phi - 0.1/(2+0.8*log10(num));
    if f >= favg
        pm = (pm1+pm2)/2-(pm1-pm2)/2*sin(pi/2*alpha);
    else
        pm = pm1;
    end
	res = pm;
end

