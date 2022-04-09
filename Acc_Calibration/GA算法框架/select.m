function res=select(popsize,individuals,sumfitness,flag)
	switch(flag)
		case 1 
			ret = Roulette_select(popsize,individuals,sumfitness);%轮盘对赌法
		case 2  
			ret = Truncation_select(popsize,individuals,sumfitness);%直接截断法
		case 3 
			ret = Modified_Roulette_select(popsize,individuals,sumfitness);%轮盘赌改进法
	end
	res = ret;
end

function ret = Roulette_select(popsize,individuals,sumfitness)
	index=[];
    sumf=individuals.fitness./sumfitness;%个体适应度比例
	for i=1:popsize   
		pick=rand;%产生一个不为0的随机概率
		while pick==0 
			pick=rand;
		end
		for j=1:popsize
			pick=pick-sumf(j);%当随机概率小于个体比例概率时，选择这个个体，然后进行下一次对赌
			if pick <= 0 
                index = [index j];
				break;
			end
		end
	end
	individuals.chrome=individuals.chrome(index,:);
	individuals.fitness=individuals.fitness(index);
	ret=individuals;%选择结束后的所有个体
end

function ret = Truncation_select(popsize,individuals)
	fitness = sort(individuals.fitness,'descend');%降排序
	individuals.fitness(4*popsize/5+1:popsize) = fitness(1:popsize/5);
	individuals.chrome(4*popsize/5+1:popsize,:) = individuals.chrome(1:popsize/5,:);
	ret=individuals;
end

function ret = Modified_Roulette_select(individuals,sumfitness)%参考文献：杨亚男《基于改进遗传算法的摄像机自标定方法》
	for j = 1: popsize/2
		pick = ceil(rand*sumfitness);
		i=1;
		while sumf < pick
			i=i+1;
			sumf = individuals.fitness(i);
			
		end	
		result(j,:) = individuals.chrome(i,:);
		index = [index,i];
		for k=1:popsize
			for m = 1:length(index)
				if k~=index(m)
					sumfitness = individuals.fitness(k);
				end
			end
		end
	end
	for k=1:popsize
			for m = 1:length(index)
				if k~=index(m)
					j=1;
					not_result(j,:) = individuals.chrome(k,:);
					j=j+1;
				end
			end
		end
	ret = result;ret2 = not_result;
end
