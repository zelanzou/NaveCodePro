function relative_erro = erro(bestchrome)
	rka =[0.99      0         0
          -0.015552 0.99012   0
          -0.01531  -0.015554 0.99024];
	 
	rbis=[0.14736
		  0.14736
          0.14736];
	g0 = 9.794731130106427;
% 	mmg=10^(-3)*g0;%µ•ŒªªªÀ„
	bx = bestchrome(1) - rbis(1);
	by = bestchrome(2) - rbis(2);
	bz = bestchrome(3) - rbis(3);
	sx = bestchrome(4) - rka(1,1);
	sy = bestchrome(5) - rka(2,2);
	sz = bestchrome(6) - rka(3,3);
	mx = bestchrome(7) - rka(2,1);
	my = bestchrome(8) - rka(3,1);
	mz = bestchrome(9) - rka(3,2);
	relative_erro = [bx,by,bz,sx,sy,sz,mx,my,mz];
    figure(3);
    plot(relative_erro);
end