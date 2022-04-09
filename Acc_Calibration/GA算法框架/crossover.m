function ret = crossover(popsize,individuals,lenchrome,LB,UB,num,flag)
	sumfitness = sum(individuals.fitness);
	fmax = max(individuals.fitness);
    favg = sumfitness/popsize;
    for i=1:popsize %是否进行交叉操作由交叉概率决定（continue控制)
		pick=rand(1,2);%产生两个[0,1]之间随机敿
		while prod(pick)==0%prod计算数组连乘，如果有一个为0则结果为0
			pick=rand(1,2);
		end
		index=ceil(pick.*popsize);%ceil向数轴右边取整，得到随机选取的两个个使,此函数中index为个体编号，也就是染色体编号
		pick=rand;%随机产生不为0的交叉概率
		while pick==0
			pick=rand;
		end
		switch(flag)
			case 1 
				pc = SGA_cross();%简单交叉算子
			case 2 
				pc =  AGA_cross(individuals,fmax,favg,index);%自适应交叉算子
			case 3 
				pc = IAGA_cross(individuals,fmax,favg,index);%线性自适应交叉算子
			case 4 
				pc = HIAGA_cross(individuals,fmax,favg,index);%异性自适应交叉算子
			case 5 
				pc = TIAGA_cross(individuals,fmax,favg,index,num);%迭代自适应交叉算子	
		end
		if pick > pc        %交叉概率决定是否交叉
			continue;           %如果随机概率大于交叉概率，那么进行下丿次循环，也就是这丿次不交叉
		end
		flag=0;                 %是否进行交叉的标志位
		while flag==0           %flag=0则进行交叉
			pick=rand;          %产生1个随机数，决定交叉位罿
			while pick==0
				pick=rand;
			end
			pos=ceil(pick*lenchrome);%随机选择交叉位置
			pick=rand;%比例因子，随机产生
			v1=individuals.chrome(index(1),pos);%被选中进行交叉的两条染色体中的第一条进行交叉的位置
			v2=individuals.chrome(index(2),pos);%被选中进行交叉的两条染色体中的第二条进行交叉的位置
			individuals.chrome(index(1),pos)=pick*v2+(1-pick)*v1;    %中间重组
			individuals.chrome(index(2),pos)=pick*v1+(1-pick)*v2;    %交叉结束
			flag1=test(lenchrome,LB,UB,individuals.chrome(index(1),:));%检验交叉后染色体的可行性
			flag2=test(lenchrome,LB,UB,individuals.chrome(index(2),:));
			if(flag1*flag2==0)
				flag=0;      %若有任意1条染色体不可行，则重新交叉
			else
				flag=1;   
			end
		end
    end
    ret=individuals.chrome;
end

function res = SGA_cross()
	pc = 0.9;
	res = pc;
end

function res = AGA_cross(individuals,fmax,favg,index)
	k1=0.9;k2=0.6; 
	f = individuals.fitness(index);
    if f(1) >= f(2)
        f_ = f(1);
    else
        f_ = f(2);
    end
    if f_ > favg
		pc = k1*(fmax-f_)/(fmax-favg);
	else
		pc = k2;
    end
	res = pc;
end

function res = IAGA_cross(individuals,fmax,favg,index)
	pc1=0.9;pc2=0.6; 
	f = individuals.fitness(index);
	if f(1) >= f(2)
        f_ = f(1);
    else
        f_ = f(2);
    end
	if f_ > favg
		pc = pc1-(pc1-pc2)*(f_-favg)/(fmax-favg);
	else
		pc = pc1;
	end
	res = pc;
end

function res = HIAGA_cross(individuals,fmax,favg,index)
	pcmax=0.9;pcmin=0.6; beta = 5.512;
	f = individuals.fitness(index);
	if f(1) >= f(2)
        f_ = f(1);
    else
        f_ = f(2);
    end
	alpha = (f_ - favg)/(fmax-favg);
	if f_ > favg
		pc = 4*(pcmax-pcmin)/(exp(2*beta*alpha)+exp(-2*beta*alpha)+2)+pcmin;
	else
		pc = pcmax;
	end
	res = pc;
end

function res  = TIAGA_cross(individuals,fmax,favg,index,num)
	phi = 0.9; pc2 = 0.6;
	f = individuals.fitness(index);
	if f(1) >= f(2)
        f_ = f(1);
    else
        f_ = f(2);
    end
	alpha = (f_ - favg)/(fmax-favg);
	pc1 = phi + 1/(2+log10(num));
    if f >= favg
        pc = (pc1+pc2)/2-(pc1-pc2)/2*sin(pi/2*alpha);
    else
        pc = pc1;
    end
	res = pc;
end