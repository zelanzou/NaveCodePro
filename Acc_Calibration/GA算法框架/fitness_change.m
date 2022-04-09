function res = fitness_change(individuals,popsize,flag)
	switch(flag)
		case 1
			fitness = linear_change(individuals,popsize);%线性变换
		case 2
			fitness = dircte_change(individuals,popsize);%直接变换
        case 3
			fitness = opposite_change(individuals,popsize);%取反变换
        case 4
			fitness = backwards_change(individuals,popsize);%取倒数变换
	end
	res = fitness;
end

function fitness = linear_change(individuals,popsize)
	sumfitness = sum(individuals.fitness);
	favg = sumfitness/popsize;
	fmin = min(individuals.fitness);
    fmax = max(individuals.fitness);
	alpha = -favg/(favg-fmin);beta = fmax*favg/(favg-fmin);
	for i = 1:popsize
		fitness(i)= alpha*individuals.fitness(i)+beta;
	end
end

function fitness = dircte_change(individuals,popsize)
	for i = 1:popsize
		fitness(i)= 50-individuals.fitness(i);
	end
end
function fitness = opposite_change(individuals,popsize)
	for i = 1:popsize
		fitness(i)= -individuals.fitness(i);
	end
end
function fitness =  backwards_change(individuals,popsize)
	for i = 1:popsize
		fitness(i)= 1-1/individuals.fitness(i);
	end
end